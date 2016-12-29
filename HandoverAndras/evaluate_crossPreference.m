clear all, close all,

addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

warning('off')
ridge = 1e-4;
userData;

for j = 1:length(user.names)
    load(['HandoverLearningOrientation_', user.names{j}, '.mat'])
    
    finalPolicyMean(j, :) = data.policyMean(end, :);
    finalPolicyCov{j} = data.policyCov{end};
end


for j =1:length(user.names)
    load(['HandoverLearningOrientation_', user.names{j}, '.mat'])
    
    hypFinal = data.hyp40(end, :); 
    sig = hypFinal(1);
    sigma2 = hypFinal(2);
    w = hypFinal(3:end); W = diag(w.^-2);
    x = data.samples;
    prefs = data.prefFeedback;
    absFeedback = data.absFeedback;
    
    Sigma = exp(-.5 * maha(x, x, W)) ;
    Sigma = Sigma + eye(size(Sigma)) * ridge;
    
    f = zeros(size(x,1), 1);
    [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
    
    % Sample a lot for policy update
    iK = eye(size(Sigma))/(Sigma);
    
    for i = 1:size(finalPolicyMean, 1)
        xsampled = finalPolicyMean(i, :);
       
        kall = exp(-.5 * maha(xsampled, x, W))
        Kxx = exp(-.5 * maha(xsampled, xsampled, W)) ;
        SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
        SigmaStar = (SigmaStar + SigmaStar')/2;
        
        ypred = kall * iK * fmap;
        
        meanFinalRew(j, i) = ypred;
        stdFinalRew(j, i) = SigmaStar.^.5;
    end
end
