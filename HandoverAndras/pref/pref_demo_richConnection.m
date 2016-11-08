clear all, close all

rew = @(x) .5*exp(- (x - 2) * 1 * (x -2)) + exp(- x * 4 * x) ;
addpath('../gp')
addpath('..')
addpath('~/svnprojects/ClassSystem/Helpers/')
addpath('../DERIVESTsuite/')
addpath('../teach-grad-hess/')

x = -2:.1:5;

for i=1:length(x)
	r(i) = rew(x(i));
end

xtest = rand(1,  40)*7 -2;
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
sig = 0.2;
w = 3.2;
figure,
for i = 10:length(xtest)
    
    ixHit = find(and(any(bsxfun(@eq, prefs(:, 1), 1:i), 2), any(bsxfun(@eq, prefs(:, 2), 1:i), 2)));
    
    
    K = exp(-.5 * maha(xtest(1:i)', xtest(1:i)', w)) ;
    K = K + eye(size(K)) * ridge;

    f = randn(i, 1)*0.0;
    fmap = nr_plgp(f, prefs(ixHit, :), K, sig);
    
%     [S, dS, ddS, beta] = grad_fmap(fmap, prefs(ixHit, :), K, sig);
    iK = eye(size(K))/(K);
%     GammaMap = ddS - iK;
%     covmap = eye(size(K))/(GammaMap + iK);
%     stds = diag(covmap).^.5;
    
    kall = exp(-.5 * maha(x(:), xtest(1:i)', w));
    ypred = kall * iK * fmap;
    stdpred = (1 - sum((kall*iK) .* kall, 2) + 2*sig^2).^.5;
   

    clf, plot(x, r)
    hold on, plot(xtest(1:i), ytest(1:i), 'r*')
    
    
    ymm = max(ytest)-min(ytest);
    fmapmm = max(fmap) - min(fmap);
    hold on, plot(x,   ymm / fmapmm* (ypred - mean(fmap)) + mean(ytest), 'k')
    hold on, plot(x,   ymm / fmapmm* (ypred+ 2*stdpred - mean(fmap)) + mean(ytest), 'k--')
    hold on, plot(x,   ymm / fmapmm* (ypred- 2*stdpred - mean(fmap)) + mean(ytest), 'k--')
%     hold on, errorbar(xtest(1:i), fmap, stds, 'r*');
    drawnow
    pause(.1)
end
