clear all, close all
% 
% addpath('../gp')
% addpath('..')
% addpath('~/svnprojects/ClassSystem/Helpers/')
% addpath('../DERIVESTsuite/')
% addpath('../teach-grad-hess/')


muContext = [0; 0];
covContext = diag([.5; .3]);

dt=.02;
tend = 2;
forceStd = 1;
numTrajs = 1;
prefsForSamples = 3;
consecutivePref = 0;
initialSamples = 20;
newSamples = 5;
epsilon = .5;
policyUpdates = 3;
ridge = 1e-4;
updateSamples = 500;
controlContexts = 30;
absFeedbackFrequency = 3; % evalutes every #... sample
siga = 3;
deterministicEvaluation = 1; % always take the most likely controller to check performance
hyperRestarts = 8;
numIterations = 5;
useCmaes = 1;
cmaesSteps = 20;

saveFileName = ['Pref', num2str(prefsForSamples), '_Consec', num2str(consecutivePref), '_AbsFreq', num2str(absFeedbackFrequency), '_InitNewSamples'...
    ,num2str(initialSamples), '_', num2str(newSamples), '_eps', num2str(round(100*epsilon)), '_siga', num2str(round(siga)), '.mat'];

try
    matlabpool 4
catch
end

xcheck = mvnrnd(muContext, covContext, controlContexts);

for iter = 1:numIterations
    
    rng('default')
    rng(iter)
    
    muK =  [ -2
        -5
        -100
        -20];
    covK = diag([1; 2; 8; 5].^2);
    
    samples = [];
    rewards = [];
    contexts = [];
    prefs = [];
    fPrior = [];
    
    %% Get initial samples
    for i = 1:initialSamples
        K = mvnrnd(muK, covK);
        x0 = mvnrnd(muContext, covContext);
        [r, xu] = simulateCartPole([x0(1) 0 x0(2) 0]', 1, dt, tend, forceStd, K);
        
        samples = [samples; K(:)'];
        contexts = [contexts; x0(:)'];
        rewards = [rewards; -log(-max(-10000, sum(r)/1e4))];
        
        if mod(length(rewards), absFeedbackFrequency) == 0
            fPrior = [fPrior; length(rewards), rewards(end) + randn(1)*siga];
        end
    end
    
    for i = 1:length(rewards)
        
        % Choose random samples
        if ~consecutivePref
            % choose random samples to compare
            randIx = randperm(length(rewards), prefsForSamples+1);
            randIx = setdiff(randIx, i);
            randIx = randIx(1:prefsForSamples);
            
        else
            % choose the previous sample to compare to
            if i > 1
                randIx = i-1;
            else
                randIx = [];
            end
        end
        
        for j = 1:length(randIx)
            
            if (rewards(i)+randn(1)*rewards(i)/100) >= (rewards(randIx(j)) + randn(1) *rewards(randIx(j))/100)
                prefs = [prefs; i randIx(j)];
            else
                prefs = [prefs; randIx(j) i];
            end
        end
    end
    A = zeros(4, 2);
    a = muK;
    
    disp(' ')
    disp(' checking initial performance ...')
    disp(' ')
    
    % check initital performance
    for i = 1:controlContexts
        
        x0 = xcheck(i, :);
        if deterministicEvaluation
            K = a + A*x0(:);
        else
            K = mvnrnd(a + A*x0(:), covK);
        end
        [r, xu] = simulateCartPole([x0(1) 0 x0(2) 0]', 1, dt, tend, forceStd, K);
        
        
        rcheck(i) = -log(-max(-10000, sum(r)/1e4));
    end
    
    RMeanSave(iter, 1) = mean(rcheck);
    RStdSave(iter, 1) = std(rcheck);
    
    
    
    %% Main loop
    for iterPolicyUpdate = 1:policyUpdates
        
        disp(' ')
        disp('New iteration, sampling... ')
        disp(' ')
        % Sample
        for i = 1:newSamples
            
            x0 = mvnrnd(muContext, covContext);
            K = mvnrnd(a + A*x0(:), covK);
            
            [r, xu] = simulateCartPole([x0(1) 0 x0(2) 0]', 1, dt, tend, forceStd, K);
            
            samples = [samples; K(:)'];
            contexts = [contexts; x0(:)'];
            rewards = [rewards; -log(-max(-10000, sum(r)/1e4))];
            
            % get absolute feedback
            if mod(length(rewards), absFeedbackFrequency) == 0
                fPrior = [fPrior; length(rewards), rewards(end) + randn(1)*siga];
            end
            
            % Choose random samples
            if ~consecutivePref
                % choose random samples to compare
                randIx = randperm(length(rewards), prefsForSamples+1);
                randIx = setdiff(randIx, length(rewards));
                randIx = randIx(1:prefsForSamples);
                
            else
                % choose the previous sample to compare to
                randIx = length(rewards)-1;
            end
            
            for j = 1:length(randIx)
                
                if (rewards(end)+randn(1)*rewards(end)/100) >= (rewards(randIx(j)) + randn(1) *rewards(randIx(j))/100)
                    prefs = [prefs; length(rewards) randIx(j)];
                else
                    prefs = [prefs; randIx(j) length(rewards)];
                end
            end
        end
        
        
        % Estimate rewards
        loghyp = log(std([samples, contexts])/2);
        
        if ~isempty(fPrior)
            dummy = quantile(fPrior(:, 2), [.25 .75]);
            sigp = (dummy(2)-dummy(1))/sqrt(2);
            sigf = std(fPrior(:, 2));
        else
            sigp = 2;
            sigf = 5;
        end
        
        loghyp = [loghyp, log([sigf sigp])]';
        loghyp_init = loghyp;
        
        disp(' ')
        disp('Optimizing reward function... ')
        disp(' ')
        
        foptsave = [];
        lhoptsave = [];
        
        objfun = @(lh) pref_loghyp_opt(lh, [samples, contexts], prefs, fPrior, ridge, siga);
        
        if useCmaes ==1
            opts.MaxFunEvals = cmaesSteps;
            parfor i = 1:hyperRestarts
                
                loghyp = log(std([samples, contexts])/2);
                
                if ~isempty(fPrior)
                    dummy = quantile(fPrior(:, 2), [.25 .75]);
                    sigp = (dummy(2)-dummy(1))/sqrt(2);
                    sigf = std(fPrior(:, 2));
                else
                    sigp = 2;
                    sigf = 5;
                end
                
                loghyp = [loghyp, log([sigf sigp])]';
                
                [~, ~, ~, ~, out] = cmaes('pref_loghyp_opt', loghyp, ones(length(loghyp), 1), opts, [samples, contexts], prefs, [], ridge, siga);
                loghyp = [out.solutions.bestever.x(:); log(siga)];
                
                
                
                foptsave(i) = out.solutions.bestever.f;
                lhoptsave(i, :) = out.solutions.bestever.x';
                
            end
        else
            
            options = optimoptions('fminunc', 'MaxFunEvals', 200);
            parfor i = 1:hyperRestarts
                
                if i > 1
                    loghypLoc = loghyp + [randn(length(loghyp)-1, 1); randn(1)*2] ;
                else
                    loghypLoc = loghyp;
                end
                
                try
                    [lhopt, fopt] = fminunc(objfun, loghypLoc, options);
                    lhoptsave(i, :) = lhopt(:)';
                    foptsave(i) = fopt;
                catch err
                    if strcmp(err.message, 'NaN in ddCdf')
                        foptsave(i) = NaN;
                        lhoptsave(i, :) = NaN*ones(1, 8);
                    else err
                        disp(err.message)
                        error('something went wrong in hyperparameter optimization')
                    end
                end
                
            end
        end
        
        % adding a standard solution
        foptsave(end+1) = objfun(loghyp_init);
        lhoptsave( end+1, :) = loghyp_init(:)';
        
        
        [~, ic] = sort(foptsave);
        
        disp('optimal values/parameters')
        [foptsave(:) lhoptsave]'
        
        
        loghyp = [lhoptsave(ic(1), :)'; log(siga)];
        hyp = exp(loghyp);
        
        w = hyp(1:6); W = diag(w.^-2);
        sigf = hyp(7);
        sigp = hyp(8);
        siga = hyp(9);
        
        
        I = zeros(length(rewards));
        Sigma = sigf^2 * exp(-.5 * maha([samples, contexts], [samples, contexts], W)) + I * ridge;
        
        activations = sum(Sigma/sigf^2, 2)/size(Sigma, 2);
        disp(['mean activation: ', num2str(mean(activations))])
        
        
        [fmap, ~, GammaMap] = nr_plgp_wPrior(zeros(length(rewards), 1), prefs, Sigma, sigp, fPrior, siga);
        
%         figure(2), clf
%         scaledRew = (rewards-min(rewards))/(max(rewards)-min(rewards));
%         scaledfmap = (fmap - min(fmap))/(max(fmap)-min(fmap));
%         [sortRew, ic] = sort(scaledRew);
%         ixStart = round(length(sortRew)*.33) + 1;
%         plot(sortRew(ixStart:end), '*'), hold on, plot(scaledfmap(ic(ixStart:end)), 'r*')
%         drawnow
        
        % Sample results to update policy
        sampUpdateContext = mvnrnd([muContext], [ covContext], updateSamples);
        for i = 1:updateSamples
            sampUpdateK(i, :) = mvnrnd(a + A*sampUpdateContext(i, :)', covK);
        end
        sampUpdate = [sampUpdateK, sampUpdateContext];
        
        kall = sigf^2 *exp(-.5 * maha([samples contexts], sampUpdate, W));
        Kxx = sigf^2 *exp(-.5 * maha(sampUpdate, sampUpdate, W)) ;
        SigmaStar = Kxx - kall' / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall;
        ypred = (kall' / Sigma) * (fmap );
        
        SigmaStar = (SigmaStar + SigmaStar')/2;
        sampUpdateRew = mvnrnd(ypred, SigmaStar);
        
        disp(' ')
        disp('Updating policy... ')
        disp(' ')
        % Update policy
        Phi = [ones(updateSamples, 1), sampUpdateContext, sampUpdateContext.^2];
        repsobjfun = @(etatheta) reps_dual_state(etatheta, epsilon, sampUpdateRew, Phi);
        options = optimset('Algorithm','active-set');
        options = optimset(options, 'GradObj','on');
        options = optimset(options, 'Display', 'off');
        etatheta = [1; zeros(5, 1)];
        etatheta = fmincon(repsobjfun, etatheta(:), diag([-1; ones(5, 1)]), [-.0001; Inf; Inf; Inf; Inf; Inf], [], [], [], [], [], options);
        
        eta = etatheta(1);
        theta = etatheta(2:end);
        
        V = Phi*theta;
        
        p = exp((sampUpdateRew' - V)/eta);
        p = p/sum(p);
        
        CC = [ones(updateSamples, 1), sampUpdateContext];
        aA = (CC' * diag(p) *CC)\CC'*diag(p)*sampUpdateK;
        
        a = aA(1, :)';
        A = aA(2:end, :)';
        
        bsxfun(@minus, sampUpdateK, a');
        dummy = sampUpdateK - bsxfun(@plus, a , A*sampUpdateContext')';
        covK = bsxfun(@times, p, dummy)'*dummy / sum(p);
        
        disp(' ')
        disp('Checking learning performance ...')
        disp(' ')
        % check performance
        for i = 1:controlContexts
            
            x0 = xcheck(i, :);
            if deterministicEvaluation
                K = a + A*x0(:);
            else
                K = mvnrnd(a + A*x0(:), covK);
            end
            [r, xu] = simulateCartPole([x0(1) 0 x0(2) 0]', 1, dt, tend, forceStd, K);
            
            
            rcheck(i) = -log(-max(-10000, sum(r)/1e4));
        end
        
        disp([' ============= UPDATE #', num2str(iterPolicyUpdate), ' ============'])
        disp(' ')
        disp([' Fixed states reward mean/std: ', num2str(mean(rcheck)), '/', num2str(std(rcheck)) ]);
        disp([' Initial Fixed states reward mean/std: ', num2str(RMeanSave(1)), '/', num2str(RStdSave(1)) ]);
        disp(' ')
        disp([' ====================================='])
        
        RMeanSave(iter, iterPolicyUpdate+1) = mean(rcheck);
        RStdSave(iter, iterPolicyUpdate+1) = std(rcheck);
        
        R.mean = RMeanSave;
        R.std = RStdSave;
        R.policyA(:, :, iter) = A;
        R.policya(:, iter) = a;
        
        
        save(saveFileName, 'R');
        
%         figure(1)
%         clf
%         subplot(2,1,1),
%         plot(RMeanSave)
%         subplot(2,1,2),
%         plot(RStdSave)
%         drawnow
        
    end
end
