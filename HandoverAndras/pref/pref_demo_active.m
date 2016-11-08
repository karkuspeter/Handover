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


for i = 2:length(xtest)
    if ytest(i) >= (ytest(i-1) + randn(1)*.05)
        prefs(i-1, :) = [i i-1];
    else
        prefs(i-1, :) = [i-1 i];
    end

    K = exp(-.5 * maha(xtest(1:i)', xtest(1:i)', w)) ;
    K = K + eye(size(K)) * ridge;

    f = randn(i, 1)*0.0;
    [fmap, ddS] = nr_plgp(f, prefs, K, sig);
    

    iK = eye(size(K))/(K);
    iGMap = eye(size(K))/(ddS - iK);
    iKG = eye(size(K))/(K + iGMap);
    
    
    if and(i > 5, mod(i, 2) == 0)
        optfun = @(x) grad_logPref(x, log([sig, w^-.5]), xtest(1:i), fmap, iKG, ridge);
        options = optimoptions('fminunc','GradObj','on');
        
        xstart = linspace(min(xtest), max(xtest), 5);
        for j = 1:length(xstart)
            optx(j) = fminunc(optfun, xstart(j), options);   
            valoptx(j) = optfun(optx(j));
        end
        
        [~, ixmax] = min(valoptx);
        
        xtest(i+1) = optx(ixmax);
        ytest(i+1) = rew(optx(ixmax));
        
        [xstart', optx', valoptx']
    end
    
    
    kall = exp(-.5 * maha(x(:), xtest(1:i)', w));
    ypred = kall * iK * fmap;
    stdpred = (1 - sum((kall*iKG) .* kall, 2)).^.5;
   

    clf, plot(x, r)
    hold on, plot(xtest(1:i), ytest(1:i), 'r*')
    
    
    ymm = max(ytest)-min(ytest);
    fmapmm = max(fmap) - min(fmap);
    subplot(2,1,1)
    hold on, plot(x,   ymm / fmapmm* (ypred - mean(fmap)) + mean(ytest), 'k')
    hold on, plot(x,   ymm / fmapmm* (ypred+ 2*stdpred - mean(fmap)) + mean(ytest), 'k--')
    hold on, plot(x,   ymm / fmapmm* (ypred- 2*stdpred - mean(fmap)) + mean(ytest), 'k--')
    hold on, plot(xtest(1:i), ymm / fmapmm* ( fmap- mean(fmap))+ mean(ytest), 'ro');
    drawnow
    subplot(2,1,2)
    
    kcurr = exp(-.5 * maha(xtest(i), xtest(1:i)', w));
    ycurr = kcurr * iK * fmap;
    stdpred = (1 + 2*sig^2 - sum((kcurr*iKG) .* kcurr, 2)).^.5;
    
    plot(x, normcdf((ypred - ycurr)/stdpred)),
    hold on, plot([xtest(i), xtest(i)], [0, 1], 'k--')
    drawnow
    pause
    
end

iKG = eye(size(K))/(K + iGMap);

% [f, df] = grad_logPref(3, log([sig, w^-.5]), xtest, fmap, iKG, ridge);
optfun = @(x) grad_logPref(x, log([sig, w^-.5]), xtest, fmap, iKG, ridge);

options = optimoptions('fminunc','GradObj','on');
optx = fminunc(optfun, 3, options);
