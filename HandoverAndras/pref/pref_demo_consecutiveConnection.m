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

ridge = 1e-4;
sig = 0.3;
w = 3.2;
figure,
for i = 2:length(xtest)
    if ytest(i) >= (ytest(i-1) + randn(1)*.05)
        prefs(i-1, :) = [i i-1];
    else
        prefs(i-1, :) = [i-1 i];
    end
    
       
    
    
    K = exp(-.5 * maha(xtest(1:i)', xtest(1:i)', w)) ;
    K = K + eye(size(K)) * ridge;

    f = randn(i, 1)*0.0;
    fmap = nr_plgp(f, prefs, K, sig);
    

    iK = eye(size(K))/(K);
    
    kall = exp(-.5 * maha(x(:), xtest(1:i)', w));
    ypred = kall * iK * fmap;
    stdpred = (1 + 2*sig^2 - sum((kall*iK) .* kall, 2)).^.5;
   

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
