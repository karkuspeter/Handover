close all, clear all
addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

warning('off')

ridge = 1e-4;

fixedActivation = 0.2;

userData;

for j =1:length(user.names)
    load(['HandoverLearningOrientation_', user.names{j}, '.mat'])
    for i= 1:size(data.hyp, 1)
        
        %% init policy
        
        absFeedback=data.absFeedback;
        absFeedback(:, 2) = (absFeedback(:, 2)-1) * 4/9 -2;
        prefs = data.prefFeedback;
        x = data.samples;
        hyp = data.hyp(i, :);
        sig = hyp(1);
        sigma2 = hyp(2);
        
        if size(data.hyp, 2) > 3
            fixedW = data.hyp(i, 3:end);
        else
            disp('have to recompute kernel weights')
            fixedW = kernelActivationTrick(x, fixedActivation);
        end
        W = diag(fixedW.^-2);
        
        Sigma = exp(-.5 * maha(x, x, W)) ;
        Sigma = Sigma + eye(size(Sigma)) * ridge;
        
        f = zeros(size(x,1), 1);
        [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
        
        iK = eye(size(Sigma))/(Sigma);
        
        %init mean
        xsampled_initPolicy = data.policyMean(1, :);
        kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
        ypred = kall * iK * fmap;
        
        
        Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
        SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
        SigmaStar = (SigmaStar + SigmaStar')/2;
        
        data.stdR_initMean(i) = sqrt(SigmaStar);
        data.meanR_initMean(i) = mean(ypred);
        
        
        %% learnt policy
        
        xsampled_initPolicy = data.policyMean(end, :);
        kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
        ypred = kall * iK * fmap;
        
        
        Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
        SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
        SigmaStar = (SigmaStar + SigmaStar')/2;
        
        data.stdR_learntMean(i) = sqrt(SigmaStar);
        data.meanR_learntMean(i) = mean(ypred);
    end
    
save(['HandoverLearningOrientation_', user.names{j}], 'data')
    
end
