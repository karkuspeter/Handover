clear all, close all

rew = @(x) 2*( .5*exp(- (x - 2) * 1 * (x -2)) + exp(- x * 4 * x)-.1174) ;
addpath('../gp')
addpath('..')
addpath('~/svnprojects/ClassSystem/Helpers/')
addpath('../DERIVESTsuite/')
addpath('../teach-grad-hess/')
addpath('../reps_demo/')

xsampled = -6:.1:9;
for i = 1:length(xsampled)
    r(i) = rew(xsampled(i));
end

Mu = 2;
Cov = 3;
epsilon = .5;

probAbsolute = 0.2;
% Get some initial preferences
absFeedback = [];
x = Mu;
for i = 2:10
    x(i) = mvnrnd(Mu, Cov, 1);
    if (rew(x(i)) + randn(1)*0.05) >= (rew(x(i-1)) + randn(1)*0.05)
        prefs(i - 1, :) = [i, i-1];
    else
        prefs(i -1 , :) = [i-1, i];
    end
    if rand(1) < probAbsolute
        absFeedback = [absFeedback; [i, rew(x(i))+randn(1)*.05]];
    end
end
i = i+1;

ridge = 1e-4;
sig = 0.1;
sigma2 = 0.2;
w = 0.55; W = w.^-2;



learning = 1;
while learning
    
    % Sample
    for j = 1:5
        x(i) = mvnrnd(Mu, Cov, 1);
        
        % Get preference
        if (rew(x(i)) + randn(1)*0.05) >= (rew(x(i-1)) + randn(1)*0.05)
            prefs(i - 1, :) = [i, i-1];
        else
            prefs(i -1 , :) = [i-1, i];
        end
        
        if rand(1) < probAbsolute
            absFeedback = [absFeedback; [i, rew(x(i))+randn(1)*.05]];
        end
        i = i + 1;
    end
    
    [f, sigpOpt, sigaOpt, quantileOpt] = pref_loghyp_gridMedianTrick(x', prefs, absFeedback, ridge, logspace(-2,1,8), logspace(-2,1,8), [0.3, 0.5, 0.7], 0);
    sig = sigpOpt;
    sigma2 = sigaOpt;
    w = medianTrick(x(:), quantileOpt); W = w.^-2;
    
    % Get latent rewards
    Sigma = exp(-.5 * maha(x(:), x(:), W)) ;
    Sigma = Sigma + eye(size(Sigma)) * ridge;
    
    f = zeros(length(x), 1);
%     [fmap, ddS] = nr_plgp(f, prefs, Sigma, sig);
    [fmap, ddS, GammaMap] = nr_plgp_wPrior(f, prefs, Sigma, sig, absFeedback, sigma2);

    % Sample a lot for policy update    
    iK = eye(size(Sigma))/(Sigma);
    
    xsampled = mvnrnd(Mu, Cov, 100);
    
    kall = exp(-.5 * maha(xsampled(:), x(:), W));
    Kxx = exp(-.5 * maha(xsampled(:), xsampled(:), W)) ;
    SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
    SigmaStar = (SigmaStar + SigmaStar')/2;

    ypred = kall * iK * fmap;
    
    Rsampled = mvnrnd(ypred, SigmaStar);
    % Update policy
    objfun = @(eta) reps_dual(eta, Rsampled, epsilon);
    options = optimset('Algorithm','active-set');
    options = optimset(options, 'GradObj','on');
    options = optimset(options, 'Display', 'off');
    eta = fmincon(objfun, 1, -1, -.01, [], [], [], [], [], options);
    p = exp(Rsampled/eta)/sum(exp(Rsampled/eta));
	
	Mu = xsampled(:)'*p(:);
	Cov = bsxfun(@minus, xsampled, Mu)'*bsxfun(@times, bsxfun(@minus, xsampled, Mu), p(:));
    
    disp([Mu, Cov])
    % Visualize
    
    xsampled = -6:.1:9;
    kall = exp(-.5 * maha(xsampled(:), x(:), W));
    Kxx = exp(-.5 * maha(xsampled(:), xsampled(:), W)) ;
    SigmaStar = ridge* eye(size(Kxx)) + Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
    SigmaStar = (SigmaStar + SigmaStar')/2;
    ypred = kall * iK * fmap;
    
    for j = 1:length(x)
        y(j) = rew(x(j));
    end
    
    stdpred = diag(SigmaStar).^.5;
    
    clf, plot(xsampled, r)
    hold on, plot(x, y, 'r*')
    if ~isempty(absFeedback)
        hold on, plot(x(absFeedback(:, 1)), absFeedback(:, 2), 'ko');
    end
    hold on, plot(xsampled, ypred, 'k');
    hold on, plot(xsampled, ypred + 2*stdpred, 'k--');
    hold on, plot(xsampled, ypred - 2*stdpred, 'k--');
    
    hold on, plot(xsampled, -exp(-.5 * (xsampled - Mu).^2 /Cov))
    pause(.25)
    
    if i > 50
        break
    end
    
end
