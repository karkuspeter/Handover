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
    
    hypFinal = data.hyp(end, :); 
    sig = hypFinal(1);
    sigma2 = hypFinal(2);
    w = hypFinal(3:end); W = diag(w.^-2);
    x = data.samples;
    prefs = data.prefFeedback;
    absFeedback = data.absFeedback ;
    absFeedback(:, 2) = (absFeedback(:, 2) - 5.5) * 4/9;
    
    Sigma = exp(-.5 * maha(x, x, W)) ;
    Sigma = Sigma + eye(size(Sigma)) * ridge;
    
    f = zeros(size(x,1), 1);
    [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
    
    % Sample a lot for policy update
    iK = eye(size(Sigma))/(Sigma);
    
    for i = 1:size(finalPolicyMean, 1)
        xsampled = finalPolicyMean(i, :);
       
        kall = exp(-.5 * maha(xsampled, x, W));
        kernelAct(j, i) = mean(kall);
        Kxx = exp(-.5 * maha(xsampled, xsampled, W)) ;
        SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
        SigmaStar = (SigmaStar + SigmaStar')/2;
        
        ypred = kall * iK * fmap;
        
        meanFinalRew(j, i) = ypred * 9/4 + 5.5;
        stdFinalRew(j, i) = (SigmaStar * (9/4)^2).^.5 ;
        
    end
    
    advRewMean(j) = mean(meanFinalRew(j, j) - meanFinalRew(j, [1:(j-1), (j+1):10]));
    advRewStd(j) = std(meanFinalRew(j, j) - meanFinalRew(j, [1:(j-1), (j+1):10]));
end
plotMatrix(meanFinalRew); ylabel('MeanFinalRew')
plotMatrix(stdFinalRew); ylabel('StdFinalRew')
plotMatrix(kernelAct); ylabel('MeanKernelAct')
