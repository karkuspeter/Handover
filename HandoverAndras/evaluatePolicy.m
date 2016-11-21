close all, clear all
addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

warning('off')

load('HandoverLearningOrientation_Andras.mat')
ridge = 1e-4;
for i = 1:size(data.policyMean, 1)
    
    absFeedback = data.absFeedback;
    absFeedback(:, 2) = (absFeedback(:, 2)-1) * 4/9 -2;
    x = data.samples;
    prefs = data.prefFeedback;
    fixedActivation = 0.2;
    fixedW = kernelActivationTrick(x, fixedActivation);
 
    loghyp_opt = data.hyp(i, :);
    
    sig = exp(loghyp_opt(1));
    sigma2 = exp(loghyp_opt(2));
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
            
    policyMean(i) = mean(ypred) * 9/4 + 2;
    policyStd(i) = std(ypred) * 9/4;
        
    xsampled_initPolicy = mvnrnd(data.policyMean(1, :), data.policyCov{1}, 10000);
    kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
    Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
    ypred = kall * iK * fmap;
    
    policyMean_orig(i) = mean(ypred) *9/4 +2;
    policyStd_orig(i) = std(ypred) *9/4;
    
    xsampled_initPolicy = data.policyMean(1, :);
    kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
    Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
    ypred = kall * iK * fmap;
    
    policyMean_origMean(i) = mean(ypred) *9/4 +2;
       
end

figure, plot(policyMean, 'k-')
hold on, plot(policyMean + 2*policyStd, 'k--')
hold on, plot(policyMean - 2*policyStd, 'k--')

hold on, plot(policyMean_orig, 'b-')
plot(policyMean_orig+2*policyStd_orig, 'b--')
plot(policyMean_orig-2*policyStd_orig, 'b--')

plot(policyMean_origMean, 'r-')


figure, plot(policyMean-policyMean_orig, 'k')
hold on, plot(policyMean-policyMean_orig + 2*(policyStd+policyStd_orig), 'k--');
hold on, plot(policyMean-policyMean_orig - 2*(policyStd+policyStd_orig), 'k--');
title('Diff orig to learned policy')

