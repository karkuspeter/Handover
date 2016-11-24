close all, clear all
addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

warning('off')

load('HandoverLearningOrientation_Peter.mat')
ridge = 1e-4;

absFeedback = data.absFeedback;
absFeedback(:, 2) = (absFeedback(:, 2)-1) * 4/9 -2;
x = data.samples;
prefs = data.prefFeedback;
fixedActivation = 0.2;
fixedW = kernelActivationTrick(x, fixedActivation);

loghyp_opt = data.hyp(end, :);

sig = exp(loghyp_opt(1));
sigma2 = exp(loghyp_opt(2));
w = fixedW; W = diag(w.^-2);

% Get latent rewards
Sigma = exp(-.5 * maha(x, x, W)) ;
kernelAct = median(mean(Sigma, 2))
Sigma = Sigma + eye(size(Sigma)) * ridge;

f = zeros(size(x,1), 1);
[fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);

% Sample a lot for policy update
iK = eye(size(Sigma))/(Sigma);


rng('default');
rng(1);

numRandSamples = 20;
randSamples = rand(numRandSamples, 4) .* (repmat(data.polparMaxLimit, numRandSamples, 1)-repmat(data.polparMinLimit, numRandSamples, 1)) + repmat(data.polparMinLimit, numRandSamples, 1);


for i = 1:size(data.policyMean, 1)
    
    
    xsampled = mvnrnd(data.policyMean(i, :), data.policyCov{i}, 10000);
    kall = exp(-.5 * maha(xsampled, x, W));
    Kxx = exp(-.5 * maha(xsampled, xsampled, W)) ;
    ypred = kall * iK * fmap;
    
    policyMean(i) = mean(ypred) * 9/4 + 5.5;
    policyStd(i) = std(ypred) * 9/4;
    
    xsampled_initPolicy = mvnrnd(data.policyMean(1, :), data.policyCov{1}, 10000);
    kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
    Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
    ypred = kall * iK * fmap;
    
    policyMean_orig(i) = mean(ypred) *9/4 + 5.5;
    policyStd_orig(i) = std(ypred) *9/4;
    
     %randSamples
     kall = exp(-.5 * maha(randSamples, x, W));
    Kxx = exp(-.5 * maha(randSamples, randSamples, W)) ;
    ypred = kall * iK * fmap;
    
    policyMean_rand(i) = mean(ypred) *9/4 +5.5;
    policyStd_rand(i) = std(ypred) *9/4;
    
    xsampled_initPolicy = data.policyMean(1, :);
    kall = exp(-.5 * maha(xsampled_initPolicy, x, W));
    Kxx = exp(-.5 * maha(xsampled_initPolicy, xsampled_initPolicy, W)) ;
    ypred = kall * iK * fmap;
    
    policyMean_origMean(i) = mean(ypred) *9/4 +5.5;
    
end

figure, plot(policyMean, 'k-')
hold on, plot(policyMean_orig, 'r-')
plot(policyMean_origMean, 'b-')
plot(policyMean_rand, 'b', 'LineWidth', 2)
legend('E[R] learned', 'E[R] init', 'R init mean', 'E[R] rand')
% hold on, plot(policyMean + 2*policyStd, 'k--')
% hold on, plot(policyMean - 2*policyStd, 'k--')
% 
% plot(policyMean_orig+2*policyStd_orig, 'r--')
% plot(policyMean_orig-2*policyStd_orig, 'r--')
% 
% plot(policyMean_rand-2*policyStd_rand, 'b--', 'LineWidth', 2)
% plot(policyMean_rand+2*policyStd_rand, 'b--', 'LineWidth', 2)

title(['Performance with full reward model'])
xlabel('policy updates')