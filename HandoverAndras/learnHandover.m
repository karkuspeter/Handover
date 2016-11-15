clear all, close all
addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

warning('off')

load('HandoverLearning_test.mat')

% Sample from policy
sample = mvnrnd(data.policyMean(end, :), data.policyCov{end}, 1);

sample = max(data.polparMinLimit, sample);
sample = min(data.polparMaxLimit, sample);

% sample = [300 1 0.8 250 20 1.2 1.2];
HandoverAndras(sample);
sample = min(max(sample, data.polparMinLimit), data.polparMaxLimit);

% get feedback
failed = input('Failed experiment? [y/N]: ', 's');
if isempty(failed)
    failed = 'n';
end


if strcmp(failed, 'y')
    % ignore sample
    disp('ignoring sample')
    data.failedExperiments = [data.failedExperiments size(data.samples,1)];
else
    % update data
    data.samples = [data.samples; sample];
    preferred = input('Is the current sample preferred over previous? [y/n/i]: ', 's');
    N = size(data.samples,1);
    if strcmp(preferred, 'y')
        if ~isempty(data.prefFeedback)
            data.prefFeedback = [data.prefFeedback; [N N-1]];
        else
            data.prefFeedback = [N N-1];
        end
    elseif strcmp(preferred, 'n')
        if ~isempty(data.prefFeedback)
            data.prefFeedback = [data.prefFeedback; [N-1 N]];
        else
            data.prefFeedback = [N-1 N];
        end
    end
    
    absfeedback = 11;
    while absfeedback > 10
        absfeedback = input('How do you rate the current sample (1-10)? [-1 to ignore]: ');
        if isempty(absfeedback)
            absfeedback = -1;
        end
    end
    if absfeedback > 0
        if ~isempty(data.absFeedback)
            data.absFeedback = [ data.absFeedback; [N, absfeedback]];
        else
            data.absFeedback = [N, absfeedback];
        end
    end
    
end

N = size(data.samples, 1);

ridge = 1e-4;
if (N >= data.initSamples) && (mod(N-data.initSamples, data.updateSamples) == 0)
    % update policy
    disp('updating policy')
    
    absFeedback = data.absFeedback;
    absFeedback(:, 2) = (absFeedback(:, 2)-1) * 4/9 -2;
    x = data.samples;
    prefs = data.prefFeedback;
    activations = 0.05:.05:.5;
    
    save_loghyp_opt = zeros(length(activations), data.numHyper);
    save_fopt = zeros(1, data.numHyper);
    
    wopt = kernelActivationTrick(data.samples, 0.05:.05:.5);
    
    
    parfor i = 1:size(wopt, 1)
        
        w = wopt(i, :);
        loghyp = log([0.5, 0.2, w]);
        % get hyperparameters
        options = optimoptions('fminunc', 'Algorithm','trust-region','GradObj','on','Hessian', 'off', 'MaxFunEvals', 100, 'TolX', 1e-3, 'TolFun', 1e-2);
        optfun = @(lh) pref_loghyp_numGrad(lh, x, prefs, absFeedback, ridge, 1);
        try
            [loghyp_opt, fopt, ~, optimOutput] = fminunc(optfun, loghyp, options);
        catch
            fopt = NaN;
            loghyp_opt = ones(1, length(loghyp))*NaN;
        end
        
        save_loghyp_opt(i, :) = loghyp_opt;
        save_fopt(i) = fopt;
        
    end
    
    for i = 1:length(save_fopt)
        w = exp(save_loghyp_opt(i, 3:end)); W = diag(w.^-2);
        Sigma = exp(-.5 * maha(x, x, W)) ;  
        median_kernel_activation(i) = median(mean(Sigma, 2));
    end
    
    
    [save_fopt_sorted, ix] = sort(save_fopt);
    
    [ save_loghyp_opt(ix, :), save_fopt_sorted', median_kernel_activation(ix)',[1:length(save_fopt)]']
    optNum = 1000;
    while optNum > length(save_fopt)
        optNum = input('Which hyper do you accept?: ');
        if isempty(optNum)
            optNum = 1;
        end
    end
    
    loghyp_opt = save_loghyp_opt(ix(optNum), :);

    sig = exp(loghyp_opt(1));
    sigma2 = exp(loghyp_opt(2));
    w = exp(loghyp_opt(3:end)); W = diag(w.^-2);
    if isempty(data.hyp)
        data.hyp = exp(loghyp_opt);
    else
        data.hyp = [data.hyp; exp(loghyp_opt)];
    end
    
    % Get latent rewards
    Sigma = exp(-.5 * maha(x, x, W)) ;
    Sigma = Sigma + eye(size(Sigma)) * ridge;
    median_kernel_activation = median(mean(Sigma, 2))
    
    plotMatrix(Sigma);
    
    f = zeros(size(x,1), 1);
    [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
    
    % Sample a lot for policy update
    iK = eye(size(Sigma))/(Sigma);
    
    xsampled = mvnrnd(data.policyMean(end, :), data.policyCov{end}, 1000);
    
    kall = exp(-.5 * maha(xsampled, x, W));
    Kxx = exp(-.5 * maha(xsampled, xsampled, W)) ;
    SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
    SigmaStar = (SigmaStar + SigmaStar')/2;
    
    ypred = kall * iK * fmap;
    
    Rsampled = mvnrnd(ypred, SigmaStar);
    
    % Update policy
    objfun = @(eta) reps_dual(eta, Rsampled', data.epsilon);
    options = optimset('Algorithm','active-set');
    options = optimset(options, 'GradObj','on');
    options = optimset(options, 'Display', 'off');
    eta = fmincon(objfun, 1, -1, -.01, [], [], [], [], [], options);
    p = exp(Rsampled/eta)/sum(exp(Rsampled/eta));
    
    Mu = xsampled'*p(:);
    Cov = bsxfun(@minus, xsampled, data.policyMean(end, :))'*bsxfun(@times, bsxfun(@minus, xsampled, data.policyMean(end, :)), p(:));
    
    data.policyMean(end+1, :) = Mu';
    data.policyCov{end+1} = Cov;
    data.policyStd(end+1, :) = diag(Cov)'.^.5;
    
end

save('HandoverLearning_test', 'data')