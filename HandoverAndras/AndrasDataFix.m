close all, clear all
addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

warning('off')

ridge = 1e-4;

fixedActivation = 0.2;
currSamples = 40;
userData;


load(['HandoverLearningOrientation_', user.names{1}, '.mat'])
for i= 1:size(data.hyp, 1)

        
        absFeedback=data.absFeedback;
        absFeedback(:, 2) = (absFeedback(:, 2)-1) * 4/9 -2;
        prefs = data.prefFeedback;
        x = data.samples;
        hyp = data.hyp(i, :);
        sig = hyp(1);
        sigma2 = hyp(2);
         fixedW = kernelActivationTrick(x, fixedActivation);
         hypExtra(i, :) = fixedW(:)';
         W = diag(fixedW.^-2);
        Sigma = exp(-.5 * maha(x, x, W)) ;
        Sigma = Sigma + eye(size(Sigma)) * ridge;
        
        f = zeros(size(x,1), 1);
        [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
        
        iK = eye(size(Sigma))/(Sigma);
        
            xsampled_initPolicy = mvnrnd(data.policyMean(1, :), data.policyCov{1}, 10000);
    kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
    Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
    ypred = kall * iK * fmap;
    
    data.meanR_init(i) = mean(ypred);
    data.stdR_init(i) = std(ypred);
        
end
    