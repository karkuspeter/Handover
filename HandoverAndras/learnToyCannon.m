clear all, close all
addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

warning('off')

load('ToyCannon.mat')

% Sample from policy
sample = mvnrnd(data.policyMean(end, :), data.policyCov{end}, 1);

% show the reward
r = simulateToyCannon_SimplePolicy2(sample, 5);

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
    
    save_loghyp_opt = zeros(data.gradRestarts, data.numHyper);
    save_fopt = zeros(1, data.numHyper);
    
    
    parfor i = 1:data.gradRestarts
        
        w = medianTrick(x, 0.5+rand(1)*2);
        loghyp = log([0.1+rand(1), 0.1+rand(1), w]);
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
    save_loghyp_opt
    save_fopt
    [~, ix] = min(save_fopt);
    loghyp_opt = save_loghyp_opt(ix(1), :)
    

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
    
    plotMatrix(Sigma);
    
    f = zeros(size(x,1), 1);
    [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);
    
    % Sample a lot for policy update
    iK = eye(size(Sigma))/(Sigma);
    
    xsampled = mvnrnd(data.policyMean(end, :), data.policyCov{end}, 100);
    
    kall = exp(-.5 * maha(xsampled, x, W));
    Kxx = exp(-.5 * maha(xsampled, xsampled, W)) ;
    SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
    SigmaStar = (SigmaStar + SigmaStar')/2;
    
    ypred = kall * iK * fmap;
    
    Rsampled = mvnrnd(ypred, SigmaStar);
    figure, plot3(xsampled(:, 1), xsampled(:, 2), ypred, '*'), grid on
    hold on, plot3(xsampled(:, 1), xsampled(:, 2), Rsampled, 'ro')
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

save('ToyCannon', 'data')