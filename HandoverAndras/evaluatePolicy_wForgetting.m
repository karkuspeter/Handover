close all, clear all
addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

warning('off')

load('HandoverLearningOrientation_Gaowei.mat')
ridge = 1e-4;

lastSamples = 40;
fixedActivation = 0.2;

rng('default');
rng(1);


for i = 1:size(data.policyMean, 1)-1
    
    currSamples = data.initSamples + (i-1) * data.updateSamples;
 
    ixHigh = currSamples;
    ixLow = max(1, currSamples-lastSamples);
    ixOk = ixLow:ixHigh;
    
    ixOkAbsFeedback = and(data.absFeedback(:, 1) >= ixLow, data.absFeedback(:, 1) <= ixHigh);
    absFeedback = data.absFeedback(ixOkAbsFeedback, :);
    absFeedback(:, 2) = (absFeedback(:, 2)-1) * 4/9 -2;
    absFeedback(:, 1) = absFeedback(:, 1) - ixLow + 1;
    x = data.samples(ixOk, :);
    prefs = data.prefFeedback;
    
    ixOkPrefs = and(and(prefs(:, 1) >= ixLow, prefs(:, 1) <=ixHigh), and(prefs(:, 2) >= ixLow, prefs(:, 2) <=ixHigh));
    prefs = prefs(ixOkPrefs, :) - ixLow +1;
    
    
    try
        fixedW = kernelActivationTrick(x, fixedActivation);
    catch
        disp('fixedW computation error, using previous one')
        fixedW = hypOpt(i-1, 3:end);
    end
    loghyp = log([0.5, 0.2]);
    % get hyperparameters
    options = optimoptions('fminunc', 'Algorithm','trust-region','GradObj','on','Hessian', 'off', 'MaxFunEvals', 1000, 'TolX', 1e-3, 'TolFun', 1e-2);
    optfun = @(lh) pref_loghyp_numGrad_fixedKernelActivation(lh, x, prefs, absFeedback, ridge, 1, fixedW);
    
    [loghyp_opt, fopt, ~, optimOutput] = fminunc(optfun, loghyp, options);
    
    sig = exp(loghyp_opt(1));
    sigma2 = min(exp(loghyp_opt(2)), 0.5);
    w = fixedW; W = diag(w.^-2);
    hypOpt(i, :) = [sig sigma2 fixedW(:)'];    
    
    % Get latent rewards
    Sigma = exp(-.5 * maha(x, x, W)) ;
    kernelAct = median(mean(Sigma, 2))
    Sigma = Sigma + eye(size(Sigma)) * ridge;
    
    f = zeros(size(x,1), 1);
    [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
    
    iK = eye(size(Sigma))/(Sigma);
    
    % learntpolicy
    xsampled = mvnrnd(data.policyMean(i, :), data.policyCov{i}, 10000);
    kall = exp(-.5 * maha(xsampled, x, W));
    Kxx = exp(-.5 * maha(xsampled, xsampled, W)) ;
    ypred = kall * iK * fmap;
    
    SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
    SigmaStar = (SigmaStar + SigmaStar')/2;
   
     Rsampled = mvnrnd(ypred, SigmaStar);
    policyStd(i) = std(Rsampled) *9/4 ;
    
    
    policyMean(i) = mean(Rsampled) * 9/4 + 5.5;
%     policyStd(i) = mean(diag(SigmaStar).^.5) * 9/4;
    
    %initpolicy
    xsampled_initPolicy = mvnrnd(data.policyMean(1, :), data.policyCov{1}, 10000);
    kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
    Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
    ypred = kall * iK * fmap;
    SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
    SigmaStar = (SigmaStar + SigmaStar')/2;
   Rsampled = mvnrnd(ypred, SigmaStar);
    policyStd_orig(i) = std(Rsampled) *9/4 ;
    
    policyMean_orig(i) = mean(Rsampled) *9/4 +5.5;
%     policyStd_orig(i) = mean(diag(SigmaStar).^.5) * 9/4;
    
    %mean Policy
    xsampled_initPolicy = data.policyMean(i, :);
    kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
    Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
    ypred = kall * iK * fmap;
    kernelActivationLearntPolicy(i) = mean(kall);
    maxKernelActivationLearntPolicy(i) = max(kall);
    medianKernelActivationLearntPolicy(i) = median(kall);
    
     
    SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
    SigmaStar = (SigmaStar + SigmaStar')/2;
   
    Rsampled = mvnrnd(ypred, SigmaStar);
    policyStd_rand(i) = std(Rsampled) *9/4 ;
    
    policyMean_rand(i) = mean(Rsampled) *9/4 +5.5;
    
%     policyStd_rand(i) = mean(diag(SigmaStar).^.5) * 9/4;
    
    %init mean
    xsampled_initPolicy = data.policyMean(1, :);
    kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
    Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
    ypred = kall * iK * fmap;
    kernelActivationInitPolicy(i) = mean(kall);
    maxKernelActivationInitPolicy(i) = max(kall);
    medianKernelActivationInitPolicy(i) = median(kall);
    
     
    SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
    SigmaStar = (SigmaStar + SigmaStar')/2;
    Rsampled = mvnrnd(ypred, SigmaStar);
    policyMean_origStd(i) = std(Rsampled) *9/4 ;
    
    
    policyMean_origMean(i) = mean(Rsampled) *9/4 +5.5;
%     policyMean_origStd(i) = mean(diag(SigmaStar).^.5) * 9/4;
end

userData;
fs = 16;
figure, plot(policyMean, 'k--')
hold on, plot(policyMean_orig, 'r-.')
plot(policyMean_rand, 'b', 'LineWidth', 2)
plot(policyMean_origMean, 'b-')
plot(14, mean(user.quizFinal(4, :)), 's', 'MarkerSize', fs)
plot(14, mean(user.quizInit(4, :)), 'o', 'MarkerSize', fs)

legend('E_{\pi}[R] learned', 'E_{\pi}[R] init', 'E[R_{\mu}] learned', 'E[R_{\mu}] init', 'Human Learned', 'Human Init' )


xlabel('policy updates')
ylabel('Reward')
title('Learning performance')
set(gca, 'FontSize', fs)
set(gca, 'box', 'off')
legend('boxoff')


figure,
hold on, plot(policyStd, 'k--') 
plot(policyStd_orig, 'r-.')
plot(policyStd_rand, 'b', 'LineWidth', 2)
plot(policyMean_origStd, 'b-');

legend('Std_{\pi}[R] learned', 'Std_{\pi}[R] init', 'Std[R_{\mu}] learned', 'Std[R_{\mu}] init' )
xlabel('policy updates')
ylabel('Standard Deviation')
title('Learning performance')
set(gca, 'FontSize', fs)
set(gca, 'box', 'off')
legend('boxoff')

% figure, plot(kernelActivationLearntPolicy, 'b', 'LineWidth', 2)
% hold on, plot(kernelActivationInitPolicy, 'r')
% plot(maxKernelActivationLearntPolicy, 'b--', 'LineWidth', 2);
% plot(maxKernelActivationInitPolicy, 'r--');
% plot(medianKernelActivationLearntPolicy, 'b.-', 'LineWidth', 1);
% plot(medianKernelActivationInitPolicy, 'r.-');
% 
% legend('kernelActivation learnt', 'kernelAcivation init', 'kernelActivation maxLearnt', 'kernelActivation MaxIniit', 'kernelActivation medianLearnt', 'kernelActivation meidanIniit')

