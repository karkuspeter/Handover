clear all, close all
addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')


warning('off')

load('HandoverLearningOrientation_test.mat')

% Sample from policy
sample = mvnrnd(data.policyMean(end, :), data.policyCov{end}, 1);

sample = max(data.polparMinLimit, sample);
sample = min(data.polparMaxLimit, sample);

HandoverAndras(sample);
sample = min(max(sample, data.polparMinLimit), data.polparMaxLimit);

% get feedback
failed = input('Failed experiment? [y/N]: ', 's');
if isempty(failed)
    failed = 'n';
end


if strcmp(failed, 'y')
    % ignore sample
    disp('ignoring sample')
    data.failedExperiments = [data.failedExperiments size(data.samples,1)];
else
    % update data
    data.samples = [data.samples; sample];
    preferred = input('Is the current sample preferred over previous? [y/n/i]: ', 's');
    N = size(data.samples,1);
    if strcmp(preferred, 'y')
        if ~isempty(data.prefFeedback)
            data.prefFeedback = [data.prefFeedback; [N N-1]];
        else
            data.prefFeedback = [N N-1];
        end
    elseif strcmp(preferred, 'n')
        if ~isempty(data.prefFeedback)
            data.prefFeedback = [data.prefFeedback; [N-1 N]];
        else
            data.prefFeedback = [N-1 N];
        end
        
    end
    
    absfeedback = 11;
    while absfeedback > 10
        absfeedback = input('How do you rate the current sample (1-10)? [-1 to ignore]: ');
        if isempty(absfeedback)
            absfeedback = -1;
        end
    end
    if absfeedback > 0
        if ~isempty(data.absFeedback)
            data.absFeedback = [ data.absFeedback; [N, absfeedback]];
        else
            data.absFeedback = [N, absfeedback];
        end
        
          
        ixLow = find(data.absFeedback(:, 2) < (absfeedback-1));
        if ~isempty(ixLow)
            prefFeedbackLow = [N*ones(length(ixLow), 1) data.absFeedback(ixLow, 1) ];
            data.prefFeedback = [data.prefFeedback; prefFeedbackLow];
        end
        
        ixHigh = find(data.absFeedback(:, 2) > (absfeedback+1));
        if ~isempty(ixHigh)
            prefFeedbackHigh = [data.absFeedback(ixHigh, 1) N*ones(length(ixHigh), 1) ];
            data.prefFeedback = [data.prefFeedback; prefFeedbackHigh];
        end
    end
    
end

N = size(data.samples, 1);

ridge = 1e-4;
if (N >= data.initSamples) && (mod(N-data.initSamples, data.updateSamples) == 0)
    % update policy
    disp('updating policy')
    
    absFeedback = data.absFeedback;
    absFeedback(:, 2) = (absFeedback(:, 2)-1) * 4/9 -2;
    x = data.samples;
    prefs = data.prefFeedback;
    fixedActivation = 0.2;
    fixedW = kernelActivationTrick(x, fixedActivation);
    
    loghyp = log([0.5, 0.2]);
    % get hyperparameters
    options = optimoptions('fminunc', 'Algorithm','trust-region','GradObj','on','Hessian', 'off', 'MaxFunEvals', 1000, 'TolX', 1e-3, 'TolFun', 1e-2);
    optfun = @(lh) pref_loghyp_numGrad_fixedKernelActivation(lh, x, prefs, absFeedback, ridge, 1, fixedW);
    
    [loghyp_opt, fopt, ~, optimOutput] = fminunc(optfun, loghyp, options);
    
    sig = exp(loghyp_opt(1));
    sigma2 = exp(loghyp_opt(2));
    w = fixedW; W = diag(w.^-2);
    if isempty(data.hyp)
        data.hyp = [exp(loghyp_opt), fixedW(:)'];
    else
        data.hyp = [data.hyp; [exp(loghyp_opt), fixedW(:)']];
    end
    
    % Get latent rewards
    Sigma = exp(-.5 * maha(x, x, W)) ;
    Sigma = Sigma + eye(size(Sigma)) * ridge;
    if isempty(data.kernelActivation)
        data.kernelActivation = median(mean(Sigma, 2));
    else
        data.kernelActivation = [data.kernelActivation, median(mean(Sigma, 2))];
    end
    
    plotMatrix(Sigma);
    
    f = zeros(size(x,1), 1);
    [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
    
    % Sample a lot for policy update
    iK = eye(size(Sigma))/(Sigma);
    
    xsampled = mvnrnd(data.policyMean(end, :), data.policyCov{end}, 1000);
    
    kall = exp(-.5 * maha(xsampled, x, W));
    Kxx = exp(-.5 * maha(xsampled, xsampled, W)) ;
    SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
    SigmaStar = (SigmaStar + SigmaStar')/2;
    
    ypred = kall * iK * fmap;
    
    Rsampled = mvnrnd(ypred, SigmaStar);
    
    % Update policy
    objfun = @(eta) reps_dual(eta, Rsampled', data.epsilon);
    options = optimset('Algorithm','active-set');
    options = optimset(options, 'GradObj','on');
    options = optimset(options, 'Display', 'off');
    eta = fmincon(objfun, 1, -1, -.01, [], [], [], [], [], options);
    p = exp(Rsampled/eta)/sum(exp(Rsampled/eta));
    
    Mu = xsampled'*p(:);
    Cov = bsxfun(@minus, xsampled, data.policyMean(end, :))'*bsxfun(@times, bsxfun(@minus, xsampled, data.policyMean(end, :)), p(:));
    
    data.policyMean(end+1, :) = Mu';
    data.policyCov{end+1} = Cov;
    data.policyStd(end+1, :) = diag(Cov)'.^.5;
    data.meanR(end+1) = mean(ypred);
    data.stdR(end+1) = std(ypred);
    
    disp('kernelActivaton')
    data.kernelActivation(end)
    figure, plot(data.meanR, 'k-')
    hold on, 
    plot(data.meanR + data.stdR*2, 'k--')
    plot(data.meanR - data.stdR*2, 'k--')
    
    % evaluate initial policy
    xsampled = mvnrnd(data.policyMean(1, :), data.policyCov{1}, 1000);
    
    kall = exp(-.5 * maha(xsampled, x, W));
    Kxx = exp(-.5 * maha(xsampled, xsampled, W)) ;
    SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
    SigmaStar = (SigmaStar + SigmaStar')/2;
    
    ypred = kall * iK * fmap;
    
    data.meanR_init(end+1) = mean(ypred);
    data.stdR_init(end+1) = std(ypred);
    
    plot(data.meanR_init, 'r-')
    plot(data.meanR_init + 2*data.stdR_init, 'r--')
    plot(data.meanR_init - 2*data.stdR_init, 'r--')
    % evaluate initial mean
    xsampled = data.policyMean(1, :);
    
    kall = exp(-.5 * maha(xsampled, x, W));
    Kxx = exp(-.5 * maha(xsampled, xsampled, W)) ;
    SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
    SigmaStar = (SigmaStar + SigmaStar')/2;
    
    ypred = kall * iK * fmap;
    data.meanR_initMean(end+1) = mean(ypred);
    plot(data.meanR_initMean, 'b-')
  
    keyboard
    
end

save('HandoverLearningOrientation_test', 'data')