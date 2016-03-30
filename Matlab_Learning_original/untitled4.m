clear all, close all

rew = @(x) .5*exp(- (x - 2) * 1 * (x -2)) + exp(- x * 4 * x) ;
addpath('../gp')
addpath('..')
addpath('~/svnprojects/ClassSystem/Helpers/')
addpath('../DERIVESTsuite/')
addpath('../teach-grad-hess/')

x = -6:.1:9;

for i=1:length(x)
	r(i) = rew(x(i));
end

xtest = rand(1,  20)*7 -2;
for i = 1:length(xtest)
	ytest(i) = rew(xtest(i));% + randn(1) * .1;
end

% Generate training data
N = 50;
for i = 1:N
   
    ixRand = randperm(length(ytest), 2);
    
    if (ytest(ixRand(1)) + randn(1)*0.05) >= (ytest(ixRand(2)) + randn(1)*0.05)
        prefs(i, :) = ixRand;
    else
        prefs(i, :) = [ixRand(2) ixRand(1)];
    end

end
ridge = 1e-4;
sig = 0.5;
w = 0.55; W = w.^-2;
figure,
for i = length(xtest)
    
    ixHit = find(and(any(bsxfun(@eq, prefs(:, 1), 1:i), 2), any(bsxfun(@eq, prefs(:, 2), 1:i), 2)));
    
    
    Sigma = exp(-.5 * maha(xtest(1:i)', xtest(1:i)', W)) ;
    Sigma = Sigma + eye(size(Sigma)) * ridge;

    f = randn(i, 1)*0.0;
    
    
    [fmap, ddS] = nr_plgp(f, prefs(ixHit, :), Sigma, sig);
    
%     [S, dS, ddS, beta] = grad_fmap(fmap, prefs(ixHit, :), Sigma, sig);
    iK = eye(size(Sigma))/(Sigma);
    GammaMap = ddS - iK;
    
    kall = exp(-.5 * maha(x(:), xtest(1:i)', W));
    Kxx = exp(-.5 * maha(x(:), x(:), W)) ;
    SigmaStar = Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
    
%     covmap = eye(size(K))/(GammaMap + iK);
%     stds = diag(covmap).^.5;
    
    kall = exp(-.5 * maha(x(:), xtest(1:i)', W));
    ypred = kall * iK * fmap;
    stdpred = (1 - sum((kall*iK) .* kall, 2)).^.5;
%     stdpred = diag(SigmaStar).^.5;
   

    clf, plot(x, r)
    hold on, plot(xtest(1:i), ytest(1:i), 'r*')
    
    
    ymm = max(ytest)-min(ytest);
    fmapmm = max(fmap) - min(fmap);
%     hold on, plot(x,   ymm / fmapmm* (ypred - mean(fmap)) + mean(ytest), 'k')
%     hold on, plot(x,   ymm / fmapmm* (ypred+ 2*stdpred - mean(fmap)) + mean(ytest), 'k--')
%     hold on, plot(x,   ymm / fmapmm* (ypred- 2*stdpred - mean(fmap)) + mean(ytest), 'k--')

    hold on, plot(x, ypred, 'k');
    hold on, plot(x, ypred + 2*stdpred, 'k--');
    hold on, plot(x, ypred - 2*stdpred, 'k--');
%     hold on, errorbar(xtest(1:i), fmap, stds, 'r*');
    drawnow
    pause(.1)
end

loghyp = log([w, 1, sig]);
optfun = @(x) optimizeQuery(x, xtest', fmap, GammaMap, loghyp, ridge);

options = optimoptions('fminunc','GradObj','on'); % indicate gradient is provided 
x = fminunc(optfun,0,options)


check = -6:.1:9;
% 
for i = 1:length(check)
    [f(i), df(i)] = optimizeQuery(check(i), xtest', fmap, GammaMap, loghyp, ridge);
    optfun = @(x) optimizeQuery(x, xtest', fmap, GammaMap, loghyp, ridge);
    ndf(i) = numericalGradient(optfun, check(i));
    
end

figure, plot3(check, df, f), hold on, plot3(check, ndf, f, '--'), grid on
xlabel('x'), ylabel('df'), zlabel('f')
