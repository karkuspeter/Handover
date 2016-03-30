clear all, close all

rew = @(x) .5*exp(- (x - 2) * 1 * (x -2)) + exp(- x * 4 * x) ;
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

Mu = 3;
Cov = 2;
epsilon = .05;

% Get some initial preferences
x = Mu;
for i = 2:10
    x(i) = mvnrnd(Mu, Cov, 1);
    if (rew(x(i)) + randn(1)*0.05) >= (rew(x(i-1)) + randn(1)*0.05)
        prefs(i - 1, :) = [i, i-1];
    else
        prefs(i -1 , :) = [i-1, i];
    end
end

ridge = 1e-4;
sig = 0.5;
w = 0.55; W = w.^-2;
learning = 1;
while learning
    
    % Sample
    for j = 1:1
        x(i) = mvnrnd(Mu, Cov, 1);
        
        % Get preference
        if (rew(x(i)) + randn(1)*0.05) >= (rew(x(i-1)) + randn(1)*0.05)
            prefs(i - 1, :) = [i, i-1];
        else
            prefs(i -1 , :) = [i-1, i];
        end
        i = i + 1;
    end
    
    % Get latent rewards
    Sigma = exp(-.5 * maha(x(:), x(:), W)) ;
    Sigma = Sigma + eye(size(Sigma)) * ridge;
    
    f = zeros(length(x), 1);
    [fmap, ddS] = nr_plgp(f, prefs, Sigma, sig);

    % Sample a lot for policy update    
    iK = eye(size(Sigma))/(Sigma);
    GammaMap = ddS - iK;
    
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
    hold on, plot(xsampled, ypred, 'k');
    hold on, plot(xsampled, ypred + 2*stdpred, 'k--');
    hold on, plot(xsampled, ypred - 2*stdpred, 'k--');
    
    hold on, plot(xsampled, -exp(-.5 * (xsampled - Mu).^2 /Cov))
    pause(.5)
    
    if i > 30
        break
    end
    
end
