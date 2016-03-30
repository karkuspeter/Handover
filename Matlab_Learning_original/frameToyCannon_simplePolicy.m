clear all, close all

addpath('../gp')
addpath('..')
% addpath('~/svnprojects/ClassSystem/Helpers/')
% addpath('../DERIVESTsuite/')
% addpath('../teach-grad-hess/')


muContext = [5];
covContext = diag([1.5]);
rangeContext = [2 10];

angleStd = 1;
numTrajs = 1;
prefsForSamples = 1;
consecutivePref = 0;
initialSamples =40;
newSamples = 5;
epsilon = .75;
policyUpdates = 30;
ridge = 1e-4;
updateSamples = 500;
controlContexts = 30;
absFeedbackFrequency = 10000; % evalutes every #... sample
siga = .1;
deterministicEvaluation = 1; % always take the most likely controller to check performance
hyperRestarts = 3;
numIterations = 15;
useCmaes = 1;
cmaesSteps = 100;

seedBias = 0;



saveFileName = ['PrefToyCannonSimple', num2str(prefsForSamples), '_Consec', num2str(consecutivePref), '_AbsFreq', num2str(absFeedbackFrequency), '_InitNewSamples'...
    ,num2str(initialSamples), '_', num2str(newSamples), '_eps', num2str(round(100*epsilon)), '_siga', num2str(round(siga)), '.mat'];

try 
	matlabpool 3
catch err
	disp(err.message)
end


rng('default'); rng(666);
% xcheck = mvnrnd(muContext, covContext, controlContexts);
xcheck = rand(controlContexts, 1) * (rangeContext(2)-rangeContext(1))+rangeContext(1);

for iter = 1:numIterations
    
    rng('default')
    rng(iter)
    
    muK =  [pi/4];
    covK = diag([.2].^2);
    
    samples = [];
    rewards = [];
    contexts = [];
    prefs = [];
    fPrior = [];
    
    %% Get initial samples
    for i = 1:initialSamples
        K = mvnrnd(muK, covK);
        %         x0 = mvnrnd(muContext, covContext);
        x0 = rand(1)*(rangeContext(2)-rangeContext(1))+rangeContext(1);
        %         [r, xu] = simulateCartPole([x0(1) 0 x0(2) 0]', 1, dt, tend, forceStd, K);
        r = simulateToyCannon_SimplePolicy(K, x0);
        
        samples = [samples; K(:)'];
        contexts = [contexts; x0(:)'];
        rewards = [rewards; r];
        
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
            
            if (rewards(i)+randn(1)*.1) >= (rewards(randIx(j)) + randn(1)*.1)
                prefs = [prefs; i randIx(j)];
            else
                prefs = [prefs; randIx(j) i];
            end
        end
    end
    A = zeros(1, 1);
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
        r = simulateToyCannon_SimplePolicy(K, x0);
        rcheck(i) = r;
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
            
            r = simulateToyCannon_SimplePolicy(K, x0);
            
            samples = [samples; K(:)'];
            contexts = [contexts; x0(:)'];
            rewards = [rewards; r];
            
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
                
                if (rewards(end)+randn(1)*.1) >= (rewards(randIx(j)) + randn(1) *.1)
                    prefs = [prefs; length(rewards) randIx(j)];
                else
                    prefs = [prefs; randIx(j) length(rewards)];
                end
            end
        end
        
        % generate artificial preferences from absolute feedback
        prefsAbs = [];
        for i = 1:size(fPrior, 1)
            ps = fPrior(i, 2) > fPrior(i+1:end, 2);
            prefsAbs = [prefsAbs; ps.* fPrior(i, 1) + (1-ps).* fPrior(i+1:end, 1), (1-ps).* fPrior(i, 1) + ps.* fPrior(i+1:end, 1)];
        end
        
        % Estimate rewards
        loghyp = log(std([samples, contexts]));
        
        if ~isempty(fPrior)
            dummy = quantile(fPrior(:, 2), [.25 .75]);
            sigp = (dummy(2)-dummy(1))/sqrt(2);
            sigf = 1;
        else
            sigp = 2;
            sigf = 3;
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
            opts.SaveVariables = 'off';
%             opts.TolFun = .1;
% profile off
% profile on
           parfor i = 1:hyperRestarts
                
               
                loghyp = log(std([samples, contexts])/2);
                
                if ~isempty(fPrior)
                    dummy = quantile(fPrior(:, 2), [.25 .75]);
                    sigp = (dummy(2)-dummy(1))/sqrt(2);
                    sigf = 3;
                else
                    sigp = 2;
                    sigf = 3;
                end
                
                
                loghyp = [loghyp, log([sigf sigp])]';
        
                
                [~, ~, ~, ~, out] = cmaes('pref_loghyp_opt', loghyp, .5*ones(length(loghyp), 1), opts, [samples, contexts], prefs, fPrior, ridge, siga);
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
%      profile viewer
%      keyboard
        
        
        % adding a standard solution
        foptsave(end+1) = objfun(loghyp_init);
        lhoptsave( end+1, :) = loghyp_init(:)';
        
        
        [~, ic] = sort(foptsave);
        
        disp('optimal values/parameters')
        [foptsave(:) lhoptsave]'
            loghyp = [lhoptsave(ic(1), :)'; log(siga)];
%         
%         ferike = [305.9924  317.2016  317.1634  316.0734  949.6917
%    -1.6047   -6.3213   -0.5325   -3.7237   -1.5795
%     2.1612    1.2783   -4.6346    1.0238   -0.6655
%    -0.5832   -3.8597    4.1181   -1.3387   -0.7185
%    -0.2059   -4.1198   -6.5078   -3.3361   -0.5934
%    -3.3181    3.5575   -2.1421    4.4047   -0.7351
%     0.8808   -1.1252    1.8351    0.4610    0.8300
%     0.7482    0.3312    1.8277    0.8206         0
%    -0.5406   -1.0809    0.4007   -0.6684    0.9188
%     6.3317    4.6167    9.2125    5.2526   -2.3026];
    
%         loghyp = ferike(2:end, 1);
    
        
%         pref_loghyp_optCV(loghyp, [samples, contexts], prefs, fPrior, ridge, siga);
        
        
        hyp = exp(loghyp);
        
        w = hyp(1:2); W = diag(w.^-2);
        sigf = hyp(3);
        sigp = hyp(4);
        siga = hyp(5);
        
        
        I = eye(length(rewards));
        Sigma = sigf^2 * exp(-.5 * maha([samples, contexts], [samples, contexts], W)) + I * ridge;
        
        activations = sum(Sigma/sigf^2, 2)/size(Sigma, 2);
        disp(['mean activation: ', num2str(mean(activations))])
        
        
        [fmap, ~, GammaMap] = nr_plgp_wPrior(zeros(length(rewards), 1), prefs, Sigma, sigp, fPrior, siga);
        
        figure(2), clf
        scaledRew = (rewards-min(rewards))/(max(rewards)-min(rewards));
        scaledfmap = (fmap - min(fmap))/(max(fmap)-min(fmap));
        
        plot(fmap, rewards, 'r*'),
        if absFeedbackFrequency <= length(fmap)
            hold on, plot(fmap(absFeedbackFrequency:absFeedbackFrequency:end), rewards(absFeedbackFrequency:absFeedbackFrequency:end), 'ko')
        end
        
        drawnow
        
        % Sample results to update policy
        %         sampUpdateContext = mvnrnd([muContext], [ covContext], updateSamples);
        sampUpdateContext = rand(updateSamples, 1) * (rangeContext(2)-rangeContext(1))+rangeContext(1);
        for i = 1:updateSamples
            sampUpdateK(i, :) = mvnrnd(a + A*sampUpdateContext(i, :)', covK);
        end
        sampUpdate = [sampUpdateK, sampUpdateContext];
        
        kall = sigf^2 *exp(-.5 * maha([samples contexts], sampUpdate, W));
        Kxx = sigf^2 *exp(-.5 * maha(sampUpdate, sampUpdate, W)) ;
        SigmaStar = Kxx - kall' / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall;
        ypred = (kall' / Sigma) * (fmap );
        
        SigmaStar = (SigmaStar + SigmaStar')/2;
	try
        	sampUpdateRew = mvnrnd(ypred, SigmaStar);
	catch
		disp('Couldn''t sample parameters, trying with diagonal covariance...')
        	sampUpdateRew = mvnrnd(ypred, diag(diag(SigmaStar)));
	end
        
        disp(' ')
        disp('Updating policy... ')
        disp(' ')
        % Update policy
        Phi = [ones(updateSamples, 1), sampUpdateContext, sampUpdateContext.^2];
        repsobjfun = @(etatheta) reps_dual_state(etatheta, epsilon, sampUpdateRew, Phi);
         options = optimoptions('fminunc','GradObj','on','Algorithm','trust-region', 'Display', 'off');
                etatheta = [log(1); randn(3, 1)];
                [etatheta, ~, flag] = fminunc(repsobjfun, etatheta(:), options);%, diag([-1; ones(3, 1)]), [-.0001;Inf; Inf; Inf], [], [], [], [], [], options);
                disp(['fminunc flag: ', num2str(flag)])
                
                eta = exp(etatheta(1));
                theta = etatheta(2:end);
        
        V = Phi*theta;
        
        p = exp((sampUpdateRew' - V)/eta);
        p = p/sum(p);
	if ~any(isnan(p))
        
        disp(['KL-div: ', num2str(sum(p.*log(p ./ (ones(length(p), 1)./ length(p)))))]);
        
        CC = [ones(updateSamples, 1), sampUpdateContext];
        aA = (CC' * diag(p) * CC)\CC'*diag(p)*sampUpdateK;
        
        a = aA(1, :)';
        A = aA(2:end, :)';
        
        bsxfun(@minus, sampUpdateK, a');
        dummy = sampUpdateK - bsxfun(@plus, a , A*sampUpdateContext')';
        covK = bsxfun(@times, p, dummy)'*dummy / sum(p);

	else 
		disp('couldn''t compute p, likely to be numerical problems due to  bad fmap, using prev. policy')
	end
        
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
             
            r = simulateToyCannon_SimplePolicy(K, x0);
    
            rcheck(i) = r;
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
	R.seed(iter) = seedBias + iter;
        
        
        save(saveFileName, 'R');
        
                figure(1)
                clf
                subplot(2,1,1),
                plot(RMeanSave)
                subplot(2,1,2),
                plot(RStdSave)
                drawnow
%         
    end
end
