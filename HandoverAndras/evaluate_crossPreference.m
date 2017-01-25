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
    absFeedback = data.absFeedback;
    absFeedback(:, 2) = (data.absFeedback(:, 2) - 5.5)/9*4;
    
    Sigma = exp(-.5 * maha(x, x, W)) ;
    Sigma = Sigma + eye(size(Sigma)) * ridge;
    
    f = zeros(size(x,1), 1);
    [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
   
    iK = eye(size(Sigma))/(Sigma);
    
    for i = 1:size(finalPolicyMean, 1)
        xsampled = finalPolicyMean(i, :);
       
        kall = exp(-.5 * maha(xsampled, x, W));
        Kxx = exp(-.5 * maha(xsampled, xsampled, W)) ;
        SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
        SigmaStar = (SigmaStar + SigmaStar')/2;
        
        ypred = kall * iK * fmap;
        
        meanFinalRew(j, i) = ypred;
        stdFinalRew(j, i) = SigmaStar.^.5;
        kernelAct(j, i) = mean(kall);
    end
    
    ix = [1:j-1, j+1:10];
    
    prefRateMean(j, :) = (meanFinalRew(j, j) - meanFinalRew(j, ix)) *9/4 ;
    prefRateStd(j, :) = (stdFinalRew(j, j)^2 + stdFinalRew(j, ix).^2).^.5 * 9/4;
end

barData(:, 1) = mean(prefRateMean , 2);
errorBarData(:, 1) = (sum(prefRateStd.^2, 2)/9^2).^.5;

barData2(:, 1) = mean(kernelAct);
errorBarData2(:,1) = 4*std(kernelAct);

figure,  h = barwitherr(4*errorBarData, barData);
xlabel('#subject')
ylabel('Diff Reward')
title('Reward Advantage')
set(gca,'box','off')
set(h(1), 'FaceColor', 'y')

figure,  barwitherr(4*errorBarData2, barData2)
ylabel('kernel activation')

corrcoef([barData, barData2])
