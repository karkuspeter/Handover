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
    
    lastSamples = 40;
    
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
        
%         fixedW = kernelActivationTrick(x, fixedActivation);
%         
%         loghyp = log([.5, 0.2]);
%         % get hyperparameters
%         options = optimoptions('fminunc', 'Algorithm','trust-region','GradObj','on','Hessian', 'off', 'MaxFunEvals', 1000, 'TolX', 1e-3, 'TolFun', 1e-2);
%         optfun = @(lh) pref_loghyp_numGrad_fixedKernelActivation(lh, x, prefs, absFeedback, ridge, 1, fixedW);
%         
%         try
%             [loghyp_opt, fopt, ~, optimOutput] = fminunc(optfun, loghyp, options);
%             sig = exp(loghyp_opt(1));
%             sigma2 = min(exp(loghyp_opt(2)), 0.5);
%         catch
%             sig = .5;
%             sigma2 = .2;
%         end
%         data.hyp40(i, :) = [sig, sigma2, w(:)'];
        hyp = data.hyp40(i, :);
        sig = hyp(1);
        sigma2 = hyp(2);
        w = hyp(3:end); W = diag(w.^-2);
        
        % Get latent rewards
        Sigma = exp(-.5 * maha(x, x, W)) ;
%         kernelAct = median(mean(Sigma, 2))
        

        Sigma = Sigma + eye(size(Sigma)) * ridge;
        
        f = zeros(size(x,1), 1);
        [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
        
        iK = eye(size(Sigma))/(Sigma);
        
        data.KLdiv(i) = KLdiv_gaussian(data.policyMean(1, :), data.policyCov{1}, data.policyMean(i, :), data.policyCov{i});
        
        xsampled_initPolicy = data.policyMean(1, :);
        kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
        data.kernelActInit(i) = mean(kall);
        data.kernelActInitMax(i) = max(kall);
        
        xsampled_currPolicy = data.policyMean(i, :);
        kall = exp(-.5 * maha(xsampled_currPolicy, x, W));
        data.kernelActCurr(i) = mean(kall);
        data.kernelActCurrMax(i) = max(kall);
        
        
%         %init mean
%         xsampled_initPolicy = data.policyMean(1, :);
%         kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
%         ypred = kall * iK * fmap;
%         kernelActivationInitPolicy(i) = mean(kall);
%         maxKernelActivationInitPolicy(i) = max(kall);
%         medianKernelActivationInitPolicy(i) = median(kall);
%         policyMean_initMean(i) = mean(ypred) *9/4 +5.5;
%         
%         
%         %% learnt policy
%         currSamples = dataSize(i);
%         
%         ixHigh = min(size(data.samples, 1), 100);
%         ixLow = max(1, ixHigh-currSamples);
%         ixOk = ixLow:ixHigh;
%         
%         ixOkAbsFeedback = and(data.absFeedback(:, 1) >= ixLow, data.absFeedback(:, 1) <= ixHigh);
%         absFeedback = data.absFeedback(ixOkAbsFeedback, :);
%         absFeedback(:, 2) = (absFeedback(:, 2)-1) * 4/9 -2;
%         absFeedback(:, 1) = absFeedback(:, 1) - ixLow + 1;
%         x = data.samples(ixOk, :);
%         prefs = data.prefFeedback;
%         
%         ixOkPrefs = and(and(prefs(:, 1) >= ixLow, prefs(:, 1) <=ixHigh), and(prefs(:, 2) >= ixLow, prefs(:, 2) <=ixHigh));
%         prefs = prefs(ixOkPrefs, :) - ixLow +1;
%         
%         fixedW = kernelActivationTrick(x, fixedActivation);
%         
%         loghyp = log([0.5, 0.2]);
%         % get hyperparameters
%         options = optimoptions('fminunc', 'Algorithm','trust-region','GradObj','on','Hessian', 'off', 'MaxFunEvals', 1000, 'TolX', 1e-3, 'TolFun', 1e-2);
%         optfun = @(lh) pref_loghyp_numGrad_fixedKernelActivation(lh, x, prefs, absFeedback, ridge, 1, fixedW);
%         
%         try
%             [loghyp_opt, fopt, ~, optimOutput] = fminunc(optfun, loghyp, options);
%             
%             sig = exp(loghyp_opt(1));
%             sigma2 = min(exp(loghyp_opt(2)), 0.5);
%         catch
%             sig = .5;
%             sigma2 = .2;
%         end
%         w = fixedW; W = diag(w.^-2);
%         
%         % Get latent rewards
%         Sigma = exp(-.5 * maha(x, x, W)) ;
%         kernelAct = median(mean(Sigma, 2))
%         Sigma = Sigma + eye(size(Sigma)) * ridge;
%         
%         f = zeros(size(x,1), 1);
%         [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
%         
%         iK = eye(size(Sigma))/(Sigma);
%         
%         %init mean
%         xsampled_initPolicy = data.policyMean(end, :);
%         kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
%         Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
%         ypred = kall * iK * fmap;
%         kernelActivationLearntPolicy(i) = mean(kall);
%         maxKernelActivationLearntPolicy(i) = max(kall);
%         medianKernelActivationLearntPolicy(i) = median(kall);
%         policyMean_learntMean(i) = mean(ypred) *9/4 +5.5;
    end
    
    data.diffPolicy = data.meanR - data.meanR_init;
    data.diffPolicyMean = data.meanR_learntMean - data.meanR_initMean;
    
    save(['HandoverLearningOrientation_', user.names{j}], 'data')
end


