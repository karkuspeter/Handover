clear all, close all

rew = @(x) 2*(.5*exp(- (x - 2) * 1 * (x -2)) + exp(- x * 4 * x)-0.1174) ;
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
w = medianTrick(xtest(:)); W = w.^-2;

for j = 1:2
    for i = 20
        
        ixHit = find(and(any(bsxfun(@eq, prefs(:, 1), 1:i), 2), any(bsxfun(@eq, prefs(:, 2), 1:i), 2)));
        
        Sigma = exp(-.5 * maha(xtest(1:i)', xtest(1:i)', W)) ;
        Sigma = Sigma + eye(size(Sigma)) * ridge;
        
        % init f
        f = randn(i, 1)*0.0;
        
        if j == 2
            ixRew = [1:1:i];
            fPrior = [ixRew(:), ytest(ixRew(:))'];
        else
            fPrior = [];
        end
        sigma2 = 0.2;
        [fmap, ddS] = nr_plgp_wPrior(f, prefs(ixHit, :), Sigma, sig, fPrior, sigma2);
        
        %check grad
        
        %         [S, dS] = grad_nr_plgp_wPrior(fmap, prefs(ixHit, :), Sigma, sig, fPrior, sigma2)
        %         optfun = @(f) grad_nr_plgp_wPrior(f, prefs(ixHit, :), Sigma, sig, fPrior, sigma2);
        %         numericalGradient(optfun, fmap)
        %         disp('checking at not MAP solution')
        %         fmap = zeros(length(f), 1);
        %         [S, dS] = grad_nr_plgp_wPrior(fmap, prefs(ixHit, :), Sigma, sig, fPrior, sigma2)
        %         optfun = @(f) grad_nr_plgp_wPrior(f, prefs(ixHit, :), Sigma, sig, fPrior, sigma2);
        %         numericalGradient(optfun, fmap)
        %         keyboard
        
        %     [S, dS, ddS, beta] = grad_fmap(fmap, prefs(ixHit, :), Sigma, sig);
        iK = eye(size(Sigma))/(Sigma);
        GammaMap = ddS - iK;
        
        kall = exp(-.5 * maha(x(:), xtest(1:i)', W));
        Kxx = exp(-.5 * maha(x(:), x(:), W)) ;
        SigmaStar = Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
        
        %     covmap = eye(size(K))/(GammaMap + iK);
        %     stds = diag(covmap).^.5;
        
        kall = exp(-.5 * maha(x(:), xtest(1:i)', W));
        ypred = kall * iK * (fmap );
        %         stdpred = (1 - sum((kall*iK) .* kall, 2)).^.5;
        stdpred = diag(SigmaStar).^.5;
        
        figure
        clf, plot(x, r)
        hold on, plot(xtest(1:i), ytest(1:i), 'r*')
        if j == 2
            hold on, plot(xtest(ixRew), ytest(ixRew), 'ko')
        end
        
        hold on, plot(x, ypred, 'k');
        hold on, plot(x, ypred + 2*stdpred, 'k--');
        hold on, plot(x, ypred - 2*stdpred, 'k--');
        
        if j == 2
            title('With prior')
        else
            title('Without prior')
        end
    end
end
