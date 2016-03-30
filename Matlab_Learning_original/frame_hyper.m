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
w = 0.4; W = w.^-2;
sigf = 1;
sigma2 = .2;

loghyp = log([w, sigf, sig, sigma2]);

for kk = 1
    for i = 20
        
        hyp = exp(loghyp);
       
        
            ixRew = 1:3:20;
            fPrior = [ixRew(:), ytest(ixRew(:))'];
      
        disp(['======== k = ', num2str(kk), ' =========='])
        disp(['current loghyp = ', num2str(loghyp(:)'), ])
        
        
        opts.MaxFunEvals = 200;
        [~, ~, ~, ~, out] = cmaes('pref_loghyp_opt', loghyp(1:3), 1*ones(3, 1), opts, xtest(:), prefs, fPrior, ridge, sigma2);
        loghyp = [out.solutions.bestever.x(:); log(sigma2)];
%         objfun = @(lh) pref_loghyp_opt(lh, xtest(:), prefs, fPrior, ridge, sigma2);
%         [f, df] = objfun(loghyp(1:3));
%         ndf = numericalGradient(objfun, loghyp(1:3));
    
        
        loghyp = fminunc(objfun, loghyp(1:end-1));
        loghyp = [loghyp(:); log(sigma2)];

        hyp = exp(loghyp);
         w = hyp(1); W = w^-2;
        sigf = hyp(2);
        sig = hyp(3);
        sigma2 = hyp(4);
        
          
        Sigma = sigf^2 * exp(-.5 * maha(xtest(1:i)', xtest(1:i)', W)) ;
        Sigma = Sigma + eye(size(Sigma)) * ridge;
        
        iK = eye(size(Sigma))/(Sigma);

        [fmap, ~, GammaMap] = nr_plgp_wPrior(zeros(length(xtest), 1), prefs, Sigma, sig, fPrior, sigma2); 
        
        disp(['new loghyp = ', num2str(loghyp(:)')])

        kall = sigf^2 *exp(-.5 * maha(x(:), xtest(1:i)', W));
        Kxx = sigf^2 *exp(-.5 * maha(x(:), x(:), W)) ;
        SigmaStar = Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
        
        kall = sigf^2 *exp(-.5 * maha(x(:), xtest(1:i)', W));
        ypred = kall * iK * (fmap );
        ypred2 = kall /(Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) *fmap;
%         stdpred = (1 - sum((kall*iK) .* kall, 2)).^.5;
            stdpred2 = diag(SigmaStar).^.5;
            stdpred = (sigf^2 -  sum((kall  * iK).* kall, 2)).^.5;
        
        figure
        clf, plot(x, r)
        hold on, plot(xtest(1:i), ytest(1:i), 'r*')
        hold on, plot(xtest(ixRew), ytest(ixRew), 'ko')
        hold on, plot(x, ypred, 'r--', 'LineWidth', 2)
     
        hold on, plot(x, ypred2, 'k');
        hold on, plot(x, ypred2 + 2*stdpred2, 'k--');
        hold on, plot(x, ypred2 - 2*stdpred2, 'k--');
        
        hold on, plot(x, ypred + 2*stdpred, 'r--');
        hold on, plot(x, ypred - 2*stdpred, 'r--');
        
            title('With prior')
        
        
       
    end
end


% 
% ridge = 1e-10;
% profile on
% [f, df] = pref_loghyp_opt(loghyp(1:3), xtest(:), prefs, fPrior, ridge, .2);
% profile viewer
% objfun = @(lh) pref_loghyp_opt(lh, xtest(:), prefs, fPrior, ridge, .2);
% numericalGradient(objfun, loghyp(1:3))

