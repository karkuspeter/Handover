clear all, close all

addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

warning('off')

load('ToyCannon5dim.mat')

ridge = 1e-4;
go = 1;
while go
    
    loghyp_opt = log(data.hyp(end, :));
    absFeedback = data.absFeedback;
    absFeedback(:, 2) = (absFeedback(:, 2)-1) * 4/9 -2; 
    x = data.samples;
    prefs = data.prefFeedback;
    
    sig = exp(loghyp_opt(1));
    sigma2 = exp(loghyp_opt(2));
    fixedActivation = 0.2;
    fixedW = kernelActivationTrick(x, fixedActivation);
    W = diag(fixedW.^-2);
    
    % Get latent rewards
    Sigma = exp(-.5 * maha(x, x, W)) ;
    Sigma = Sigma + eye(size(Sigma)) * ridge;
    median_kernel_activation = median(mean(Sigma, 2));
    
%     plotMatrix(Sigma);
    
    f = zeros(size(x,1), 1);
    [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
    
    % Sample a lot for policy update
    iK = eye(size(Sigma))/(Sigma);
    
    xsampled = mvnrnd(data.policyMean(end, :), data.policyCov{end}, 5000);
    
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
    
    if all(data.policyStd(end, :) < [.01, .01, .01, .01, .01])
        go = 0;
    end
    
    try
        data.Rsave = [data.Rsave, mean(ypred)];
        data.RsaveStd = [data.RsaveStd, std(ypred)];
    catch
        data.Rsave = [mean(ypred)];
        data.RsaveStd = std(ypred);
    end
           
    
    [data.policyMean(end, :); data.policyStd(end, :)]
    clf, 
    subplot(3,1,1)
    plot(log(abs(data.policyMean)))
    ylabel('log policymean')
    subplot(3,1,2)
    plot(log(data.policyStd));
    ylabel('log policyStd')
    subplot(3,2,5)
    plot(data.Rsave)
    ylabel('log meanRew')
    subplot(3,2,6)
    plot(data.RsaveStd)
    ylabel('log stdRew')
    
    drawnow
end