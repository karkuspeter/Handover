clear all, close all
addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')



load('SimpleLearning.mat')
rew = @(x) 2*( .5*exp(- (x - 2) * 1 * (x -2)) + exp(- x * 4 * x)-.1174) ;


% Sample from policy
sample = mvnrnd(data.policyMean(end, :), data.policyCov{end}, 1);

% show the reward
disp(['Reward = ', num2str(rew(sample) + mvnrnd(0, 1)*.1)])

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
    
    % get hyperparameters
    [f, sigpOpt, sigaOpt, quantileOpt] = pref_loghyp_gridMedianTrick(x, prefs, absFeedback, ridge, logspace(-2,1,8), logspace(-2,1,8), [0.3, 0.5, 0.7], 1);
    sig = sigpOpt;
    sigma2 = sigaOpt;
    w = medianTrick(x, quantileOpt); W = w.^-2;
    if isempty(data.hyp)
        data.hyp = [sigpOpt, sigaOpt, quantileOpt];
    else
        data.hyp = [data.hyp; [sigpOpt, sigaOpt, quantileOpt]];
    end
    
    % Get latent rewards
    Sigma = exp(-.5 * maha(x, x, W)) ;
    Sigma = Sigma + eye(size(Sigma)) * ridge;
    
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
    
end

save('SimpleLearning', 'data')