close all, clear all
addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

warning('off')

ridge = 1e-4;

fixedActivation = 0.2;
dataSize = [30 50 70 100];

userData;

for j =length(user.names):length(user.names)
    load(['HandoverLearningOrientation_', user.names{j}, '.mat'])
    for i = 1:length(dataSize)
        %% init policy
        currSamples = dataSize(i);
        
        ixHigh = min(currSamples, size(data.samples, 1));
        ixLow = 1; 
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

        loghyp = log([.5, 0.2]);
        % get hyperparameters
        options = optimoptions('fminunc', 'Algorithm','trust-region','GradObj','on','Hessian', 'off', 'MaxFunEvals', 1000, 'TolX', 1e-3, 'TolFun', 1e-2);
        optfun = @(lh) pref_loghyp_numGrad_fixedKernelActivation(lh, x, prefs, absFeedback, ridge, 1, fixedW);
        
        try
        [loghyp_opt, fopt, ~, optimOutput] = fminunc(optfun, loghyp, options);
        
        sig = exp(loghyp_opt(1));
        sigma2 = min(exp(loghyp_opt(2)), 0.5);
        catch 
            sig = .5;
            sigma2 = .2;
        end
        w = fixedW; W = diag(w.^-2);
        
        % Get latent rewards
        Sigma = exp(-.5 * maha(x, x, W)) ;
        kernelAct = median(mean(Sigma, 2))
        Sigma = Sigma + eye(size(Sigma)) * ridge;
        
        f = zeros(size(x,1), 1);
        [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
        
        iK = eye(size(Sigma))/(Sigma);
        
        %init mean
        xsampled_initPolicy = data.policyMean(1, :);
        kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
        ypred = kall * iK * fmap;
        kernelActivationInitPolicy(i) = mean(kall);
        maxKernelActivationInitPolicy(i) = max(kall);
        medianKernelActivationInitPolicy(i) = median(kall);
        policyMean_initMean(i) = mean(ypred) *9/4 +5.5;
        %% learnt policy
        currSamples = dataSize(i);
        
        ixHigh = min(size(data.samples, 1), 100);
        ixLow = max(1, ixHigh-currSamples); 
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
        
      try
        [loghyp_opt, fopt, ~, optimOutput] = fminunc(optfun, loghyp, options);
        
        sig = exp(loghyp_opt(1));
        sigma2 = min(exp(loghyp_opt(2)), 0.5);
        catch 
            sig = .5;
            sigma2 = .2;
        end
        w = fixedW; W = diag(w.^-2);
                
        % Get latent rewards
        Sigma = exp(-.5 * maha(x, x, W)) ;
        kernelAct = median(mean(Sigma, 2))
        Sigma = Sigma + eye(size(Sigma)) * ridge;
        
        f = zeros(size(x,1), 1);
        [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
        
        iK = eye(size(Sigma))/(Sigma);
        
        %init mean
        xsampled_initPolicy = data.policyMean(end, :);
        kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
        Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
        ypred = kall * iK * fmap;
        kernelActivationLearntPolicy(i) = mean(kall);
        maxKernelActivationLearntPolicy(i) = max(kall);
        medianKernelActivationLearntPolicy(i) = median(kall);
        policyMean_learntMean(i) = mean(ypred) *9/4 +5.5;
        
    end
    
    figure, plot(dataSize, policyMean_learntMean, 'b')
    hold on, plot(dataSize, policyMean_initMean, 'r')
    plot(dataSize, repmat(mean(user.initScores(j, :)), 1, length(dataSize)), 'r--')
    plot(dataSize, repmat(mean(user.finalScores(j, :)), 1, length(dataSize)), 'b--')
    
    legend('learnt', 'init')
    title(user.names{j})
    ylabel('E[R]')
    xlabel('numLastSamples')
    
    saveas(gca, ['evaluations/numLastSamples_mean_', user.names{j},'.fig'])
    
    figure, plot(dataSize, kernelActivationLearntPolicy, 'b', 'LineWidth', 2)
    hold on, plot(dataSize, kernelActivationInitPolicy, 'r')
    plot(dataSize, maxKernelActivationLearntPolicy, 'b--', 'LineWidth', 2);
    plot(dataSize, maxKernelActivationInitPolicy, 'r--');
    plot(dataSize,medianKernelActivationLearntPolicy, 'b.-', 'LineWidth', 1);
    plot(dataSize, medianKernelActivationInitPolicy, 'r.-');
    xlabel('numLastSamples')
    ylabel('kernelActivation')
    title(user.names{j})
    legend('kernelActivation learnt', 'kernelAcivation init', 'kernelActivation maxLearnt', 'kernelActivation MaxIniit', 'kernelActivation medianLearnt', 'kernelActivation meidanIniit')
    
    saveas(gca, ['evaluations/numLastSamples_kernel_', user.names{j},'.fig'])
    
    
end
