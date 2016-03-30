clear all, close all

rng('default');


A = zeros(7, 1);
a = [2.75; 275; 450; 275; -2.75; 200; 550];
minmaxRanges = [0.5 5;
    50 500;
    100 800;
    50 500;
    -5 -.5
    100 300
    300 800];
covK = diag([1.12^2, 112^2, 125^2, 112^2, 1.12^2, 50^2, 125^2]);

% aAC=    [    2.7467    0.0959    1.0519
%   390.4763  -25.4907   77.2451
%   409.4036  -22.0363  109.1298
%   290.0621   -7.0848   65.3234
%    -2.8337    0.0082    0.8883
%   191.7182    6.5552   42.9994
%   388.3307   23.7705   84.9906];
%
% a = aAC(:, 1);
% A = aAC(:, 2);
% covK = diag(aAC(:, 3).^2);

steps = 50;
epsilon = .5;

loghyp = [   -0.1229
    5.7932
    8.6764
    5.4606
    1.7476
    5.2513
    5.6796
    2.9351
    1.5827
    0.4467; log(1)];

updateSamples = 500;
ridge = 1e-4;

hyp = exp(loghyp);

w = hyp(1:end-3); W = diag(w.^-2);
sigf = hyp(end-2);
sigp = hyp(end-1);
siga = hyp(end);

rng(666);
% get samples
[fPrior, prefs, context, samples] = humanFeedbackRobotHandover();
N = size(samples, 1);
I = eye(N);
    Sigma = sigf^2 * exp(-.5 * maha([samples, context], [samples, context], W)) + I * ridge;
[fmap, ddS, GammaMap] = nr_plgp_wPrior(zeros(size(samples, 1), 1), prefs, Sigma, sigp, fPrior, siga);

ix1 = find(context == 1);
ix2 = find(context == 2);
ix3 = find(context == 3);
ix4 = find(context == 4);

ixSeen = [];

for iter = 1:20
    rng('default')
    rng(iter)
    
    rewMean = []; rewStd = [];
    aAs = [];
    A = zeros(7, 1);
a = [2.75; 275; 450; 275; -2.75; 200; 550];
minmaxRanges = [0.5 5;
    50 500;
    100 800;
    50 500;
    -5 -.5
    100 300
    300 800];
covK = diag([1.12^2, 112^2, 125^2, 112^2, 1.12^2, 50^2, 125^2]);

for e = 1:steps
    
    % "get" new samples
    
    % compute likelihood for all samples
    
    dummy = bsxfun(@plus, a, A* context')'-samples;
    
    ll = exp(-sum((dummy / covK).*dummy, 2));%+.001;
    
    if e == 1
        updateSamples = 40;
    else
        updateSamples = 10;
    end
    
    ixSeenNew = [];
    for i = 1:updateSamples
        contextLoc = randperm(4, 1);
        ixCont = find(context == contextLoc);
        newIx = find([0; cumsum(ll(ixCont))/sum(ll(ixCont))] <= rand(1));
        ixSeen = [ixSeen; ixCont(newIx(end))];
        ixSeenNew(i) = ixCont(newIx(end));
    end
    
    
    % get random contexts and add the most likely ot the seen indices
    
    
    
    Nseen = length(ixSeen);
    
    
    
    % get reward model
    
    
    I = eye(Nseen);
    Sigma = sigf^2 * exp(-.5 * maha([samples(ixSeen, :), context(ixSeen, :)], [samples(ixSeen, :), context(ixSeen, :)], W)) + I * ridge;
    
    activations = sum(Sigma/sigf^2, 2)/size(Sigma, 2);
    disp(['mean activation: ', num2str(mean(activations))])
    
    
    %     [fmap, ddS, GammaMap] = nr_plgp_wPrior(zeros(size(samples, 1), 1), prefs, Sigma, sigp, fPrior, siga);
    
    
    
    
    % evaluate artifical samples
    
    for i = 1:500
        sampUpdateContext(i) = randperm(4, 1);
        sampUpdateK(i, :) = max(min(mvnrnd(a + A*sampUpdateContext(i), covK), minmaxRanges(:, 2)'), minmaxRanges(:, 1)');
    end
    sampUpdateContext =sampUpdateContext(:);
    sampUpdate = [sampUpdateK, sampUpdateContext];
    
    
    
    kall = sigf^2 *exp(-.5 * maha([samples(ixSeen, :) context(ixSeen, :)], sampUpdate, W));
    Kxx = sigf^2 *exp(-.5 * maha(sampUpdate, sampUpdate, W)) ;
    SigmaStar = Kxx - kall' / (Sigma + eye(Nseen)/(GammaMap(ixSeen, ixSeen) + ridge*eye(Nseen))) * kall;
    ypred = (kall' / Sigma) * (fmap(ixSeen, :) );
    
    SigmaStar = (SigmaStar + SigmaStar')/2;
    sampUpdateRew = ypred;%mvnrnd(ypred, SigmaStar);
    
    % update
      rewMean(e) = mean(fmap(ixSeenNew));
    rewStd(e) = std(fmap(ixSeenNew));
    aAs(:, :, e) = [a, A];
    

    ixYes = [];
    for i = 1:length(ixSeenNew)
        peti = find(fPrior(:, 1) == ixSeenNew(i));
        if ~isempty(peti)
            ixYes = [ixYes; peti];
        end
    end
    S.rewMeanAbs{iter} = mean(fPrior(ixYes, 2));
    S.rewMean{iter} = rewMean;
        S.rewStd{iter} = rewStd;
        S.aA{iter} = aAs;
        save Shandover.mat S
    if sum(ll) < 1e-3
        break
    end
        
    
    
    disp(' ')
    disp('Updating policy... ')
    disp(' ')
    % Update policy
    Phi = [ones(500, 1), sampUpdateContext, sampUpdateContext.^2];
    repsobjfun = @(etatheta) reps_dual_state(etatheta, epsilon, sampUpdateRew, Phi);
    options = optimoptions('fminunc','GradObj','on','Algorithm','trust-region', 'Display', 'off');
    etatheta = [log(1); randn(3, 1)];
    [etatheta, ~, flag] = fminunc(repsobjfun, etatheta(:), options);%, diag([-1; ones(3, 1)]), [-.0001;Inf; Inf; Inf], [], [], [], [], [], options);
    disp(['fminunc flag: ', num2str(flag)])
    
    eta = exp(etatheta(1));
    theta = etatheta(2:end);
    
    V = Phi*theta;
    
    p = exp((sampUpdateRew - V)/eta);
    p = p/sum(p);
    
    if any(isnan(p))
        break
    end
    
    disp(['KL-div: ', num2str(sum(p.*log(p ./ (ones(length(p), 1)./ length(p)))))]);
    
    CC = [ones(500, 1), sampUpdateContext];
    aA = (CC' * diag(p) * CC)\CC'*diag(p)*sampUpdateK;
    
    a = aA(1, :)';
    A = aA(2:end, :)';
    
    bsxfun(@minus, sampUpdateK, a');
    dummy = sampUpdateK - bsxfun(@plus, a , A*sampUpdateContext')';
    covK = bsxfun(@times, p, dummy)'*dummy / sum(p);
    
    disp(' ')
    disp('Checking learning performance ...')
    disp(' ')
    disp([mean(fmap(ixSeenNew)), std(fmap(ixSeenNew))])
    disp([a, A])
  
    
    
end
end