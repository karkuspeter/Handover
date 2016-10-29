clear all, close all,


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

% loghyp = [1.8162
%     6.2621
%     5.7832
%     4.4546
%     4.0078
%     6.9090
%     7.8238
%     2.6737
%     1.9727
%     1.0695
%     log(1)];

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

for e = 1:steps
    
    % get samples
    [fPrior, prefs, context, samples] = humanFeedbackRobotHandover();
    
    % get reward model
        hyp = exp(loghyp);
    
    w = hyp(1:end-3); W = diag(w.^-2);
    sigf = hyp(end-2);
    sigp = hyp(end-1);
    siga = hyp(end);
    
    I = eye(size(samples, 1));
    Sigma = sigf^2 * exp(-.5 * maha([samples, context], [samples, context], W)) + I * ridge;
    
    activations = sum(Sigma/sigf^2, 2)/size(Sigma, 2);
    disp(['mean activation: ', num2str(mean(activations))])
    
    
    [fmap, ddS, GammaMap] = nr_plgp_wPrior(zeros(size(samples, 1), 1), prefs, Sigma, sigp, fPrior, siga);
    
    

    % evaluate artifical samples
     
        for i = 1:updateSamples
            sampUpdateContext(i) = randperm(4, 1);
            sampUpdateK(i, :) = max(min(mvnrnd(a + A*sampUpdateContext(i), covK), minmaxRanges(:, 2)'), minmaxRanges(:, 1)');
        end
        sampUpdateContext =sampUpdateContext(:); 
        sampUpdate = [sampUpdateK, sampUpdateContext];
        
        
        
        kall = sigf^2 *exp(-.5 * maha([samples context], sampUpdate, W));
        Kxx = sigf^2 *exp(-.5 * maha(sampUpdate, sampUpdate, W)) ;
        SigmaStar = Kxx - kall' / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall;
        ypred = (kall' / Sigma) * (fmap );
        
        SigmaStar = (SigmaStar + SigmaStar')/2;
        sampUpdateRew = mvnrnd(ypred, SigmaStar);
    
    % update
    
    
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
    
    disp(['KL-div: ', num2str(sum(p.*log(p ./ (ones(length(p), 1)./ length(p)))))]);
    
    CC = [ones(updateSamples, 1), sampUpdateContext];
    aA = (CC' * diag(p) * CC)\CC'*diag(p)*sampUpdateK;
    
    a = aA(1, :)';
    A = aA(2:end, :)';
    
    bsxfun(@minus, sampUpdateK, a');
    dummy = sampUpdateK - bsxfun(@plus, a , A*sampUpdateContext')';
    covK = bsxfun(@times, p, dummy)'*dummy / sum(p);
    
    disp(' ')
    disp('Checking learning performance ...')
    disp(' ')
    disp([mean(ypred), std(ypred)])
    disp([a, A])
    
    aAsave(:, :, e) = [a A];
    Rsave(e) = mean(ypred);
    
    keyboard
   
end