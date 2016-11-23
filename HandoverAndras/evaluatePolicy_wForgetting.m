close all, clear all
addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

warning('off')

load('HandoverLearningOrientation_Ziquan.mat')
ridge = 1e-4;

lastSamples = 49;
fixedActivation = 0.2;

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
    
    
    fixedW = kernelActivationTrick(x, fixedActivation);
    loghyp = log([0.5, 0.2]);
    % get hyperparameters
    options = optimoptions('fminunc', 'Algorithm','trust-region','GradObj','on','Hessian', 'off', 'MaxFunEvals', 1000, 'TolX', 1e-3, 'TolFun', 1e-2);
    optfun = @(lh) pref_loghyp_numGrad_fixedKernelActivation(lh, x, prefs, absFeedback, ridge, 1, fixedW);
    
    [loghyp_opt, fopt, ~, optimOutput] = fminunc(optfun, loghyp, options);
    
    sig = exp(loghyp_opt(1));
    sigma2 = exp(loghyp_opt(2));
    w = fixedW; W = diag(w.^-2);
    hypOpt(i, :) = [sig sigma2 fixedW(:)'];    
    
    % Get latent rewards
    Sigma = exp(-.5 * maha(x, x, W)) ;
    kernelAct = median(mean(Sigma, 2))
    Sigma = Sigma + eye(size(Sigma)) * ridge;
    
    f = zeros(size(x,1), 1);
    [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
    
    % Sample a lot for policy update
    iK = eye(size(Sigma))/(Sigma);
    
    xsampled = mvnrnd(data.policyMean(i, :), data.policyCov{i}, 10000);
    kall = exp(-.5 * maha(xsampled, x, W));
    Kxx = exp(-.5 * maha(xsampled, xsampled, W)) ;
    ypred = kall * iK * fmap;
    
    policyMean(i) = mean(ypred) * 9/4 + 5.5;
    policyStd(i) = std(ypred) * 9/4;
    
    xsampled_initPolicy = mvnrnd(data.policyMean(1, :), data.policyCov{1}, 10000);
    kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
    Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
    ypred = kall * iK * fmap;
    
    policyMean_orig(i) = mean(ypred) *9/4 +5.5;
    policyStd_orig(i) = std(ypred) *9/4;
    
    xsampled_initPolicy = data.policyMean(1, :);
    kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
    Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
    ypred = kall * iK * fmap;
    
    policyMean_origMean(i) = mean(ypred) *9/4 +5.5;
    
end

figure, plot(policyMean, 'k-')
hold on, plot(policyMean_orig, 'r-')
plot(policyMean_origMean, 'b-')
legend('E[R] learned', 'E[R] init', 'R init mean')
hold on, plot(policyMean + 2*policyStd, 'k--')
hold on, plot(policyMean - 2*policyStd, 'k--')

plot(policyMean_orig+2*policyStd_orig, 'r--')
plot(policyMean_orig-2*policyStd_orig, 'r--')

title(['Performance with last ', num2str(lastSamples) , ' samples'])
xlabel('policy updates')

% Evaluate all performance with last hyp
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
    
    sig = hypOpt(end, 1);
    sigma2 = hypOpt(end, 2);
%     w = hypOpt(end, 3:end); W = diag(w.^-2);
    fixedW = kernelActivationTrick(x, fixedActivation);
    w = fixedW; W = diag(w.^-2);
    
    % Get latent rewards
    Sigma = exp(-.5 * maha(x, x, W)) ;
    kernelAct = median(mean(Sigma, 2))
    Sigma = Sigma + eye(size(Sigma)) * ridge;
    
    f = zeros(size(x,1), 1);
    [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
    
    % Sample a lot for policy update
    iK = eye(size(Sigma))/(Sigma);
    
    
    xsampled = mvnrnd(data.policyMean(i, :), data.policyCov{i}, 10000);
    kall = exp(-.5 * maha(xsampled, x, W));
    Kxx = exp(-.5 * maha(xsampled, xsampled, W)) ;
    ypred = kall * iK * fmap;
    
    policyMean(i) = mean(ypred) * 9/4 + 5.5;
    policyStd(i) = std(ypred) * 9/4;
    
    xsampled_initPolicy = mvnrnd(data.policyMean(1, :), data.policyCov{1}, 10000);
    kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
    Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
    ypred = kall * iK * fmap;
    
    policyMean_orig(i) = mean(ypred) *9/4 +5.5;
    policyStd_orig(i) = std(ypred) *9/4;
    
    xsampled_initPolicy = data.policyMean(1, :);
    kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
    Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
    ypred = kall * iK * fmap;
    
    policyMean_origMean(i) = mean(ypred) *9/4 +5.5;
    
end

figure, plot(policyMean, 'k-')
hold on, plot(policyMean_orig, 'r-')
plot(policyMean_origMean, 'b-')
legend('E[R] learned', 'E[R] init', 'R init mean')
hold on, plot(policyMean + 2*policyStd, 'k--')
hold on, plot(policyMean - 2*policyStd, 'k--')

plot(policyMean_orig+2*policyStd_orig, 'r--')
plot(policyMean_orig-2*policyStd_orig, 'r--')

title(['Performance with last ', num2str(lastSamples) , ' samples (fixed last hyp)'])
xlabel('policy updates')

