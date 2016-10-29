clear all, close all

rew = @(x) 2*(.5*exp(- (x - 2) * 1 * (x -2)) + exp(- x * 4 * x) -0.1174);
addpath('../gp')
addpath('..')
addpath('~/svnprojects/ClassSystem/Helpers/')
addpath('../DERIVESTsuite/')
addpath('../teach-grad-hess/')

x = -6:.1:9;

for i=1:length(x)
    r(i) = rew(x(i));
end

sampleNum = 50;

xtest = rand(1,  sampleNum)*7 -2;
for i = 1:length(xtest)
    ytest(i) = rew(xtest(i)) + randn(1) * .1;
end

% Generate training data
N =200;
for i = 1:N
    
    ixRand = randperm(length(ytest), 2);
    
    if (ytest(ixRand(1)) ) >= (ytest(ixRand(2)))
        prefs(i, :) = ixRand;
    else
        prefs(i, :) = [ixRand(2) ixRand(1)];
    end
    
end
ridge = 1e-4;


for kk = 1
    for i = sampleNum
        
      
            ixRew = 1:2:i;
            fPrior = [ixRew(:), ytest(ixRew(:))'];
            
%           fPrior = [];
%           ixRew = [];
            
        
    optfun = @(logsig) pref_loghyp_MedianTrick(logsig, xtest', prefs, fPrior, ridge);
    %%
%     [xx, yy]=meshgrid(0.1:.1:5, 0.1: .1:5);
%     
%     xxC = xx(:);
%     yyC = yy(:);
%     for jani = 1:length(xxC)
%         zzC(jani) = optfun(log([xxC(jani), yyC(jani)]));
%     end
%     figure, contourf(reshape(xxC, size(xx,1), size(xx,2)), reshape(yyC, size(yy, 1), size(yy, 2)), reshape(zzC, size(xx,1), size(xx,2)));
%     
%     
%     keyboard
    
%%
%     tic
%      [logsigopt, fopt] = fminunc(optfun, log([1, 1, 0.5]));
%     fminuncTime = toc;
%     
%             w = medianTrick(xtest', exp(logsigopt(3))); W = w.^-2;
%         sigf = 1;
%         sig = exp(logsigopt(1));
%         sigma2 = exp(logsigopt(2));
%         
%           
%         Sigma = sigf^2 * exp(-.5 * maha(xtest(1:i)', xtest(1:i)', W)) ;
%         Sigma = Sigma + eye(size(Sigma)) * ridge;
%         
%         iK = eye(size(Sigma))/(Sigma);
% 
%         [fmap, ~, GammaMap] = nr_plgp_wPrior(zeros(length(xtest), 1), prefs, Sigma, sig, fPrior, sigma2); 
%         
% 
%         kall = sigf^2 *exp(-.5 * maha(x(:), xtest(1:i)', W));
%         Kxx = sigf^2 *exp(-.5 * maha(x(:), x(:), W)) ;
%         SigmaStar = Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
%         
%                 
%         ypred2 = kall /(Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) *fmap;
% %         stdpred = (1 - sum((kall*iK) .* kall, 2)).^.5;
%             stdpred2 = diag(SigmaStar).^.5;
%          
%         
%         figure
%         clf, plot(x, r)
%         hold on, plot(xtest(1:i), ytest(1:i), 'r*')
%         hold on, plot(xtest(ixRew), ytest(ixRew), 'ko')
%         
%      
%         hold on, plot(x, ypred2, 'k');
%         hold on, plot(x, ypred2 + 2*stdpred2, 'k--');
%         hold on, plot(x, ypred2 - 2*stdpred2, 'k--');
%         
%             title('Fminunc')
        
        
    
    
%%
    resolutionSigmaP = logspace(-2,1,10);
    resolutionSigmaA= logspace(-2,1,10);
    resolutionQuantile = 0.3:0.1:0.7;
    
%     resolutionSigmaP = 10.1;
%     resolutionSigmaA = 0.1;
    resolutionQuantile = 0.5;
    
    tic
     [f, sigpOpt, sigaOpt, quantileOpt] = pref_loghyp_gridMedianTrick(xtest', prefs, fPrior, ridge, resolutionSigmaP, resolutionSigmaA, resolutionQuantile, 1);
     gridTime = toc;
%     figure, contourf(f, 30); xlabel('sigmaA'), ylabel('sigmaP'), title(['Opt sigP/sigA = ', num2str(sigpOpt), '/', num2str(sigaOpt)]);
    
    
        w = medianTrick(xtest',quantileOpt); W = w.^-2;
        sigf = 1;
        sig = sigpOpt;
        sigma2 = sigaOpt;
        
          
        Sigma = sigf^2 * exp(-.5 * maha(xtest(1:i)', xtest(1:i)', W)) ;
        Sigma = Sigma + eye(size(Sigma)) * ridge;
        
        iK = eye(size(Sigma))/(Sigma);

        [fmap, ~, GammaMap] = nr_plgp_wPrior(zeros(length(xtest), 1), prefs, Sigma, sig, fPrior, sigma2); 
        

        kall = sigf^2 *exp(-.5 * maha(x(:), xtest(1:i)', W));
        Kxx = sigf^2 *exp(-.5 * maha(x(:), x(:), W)) ;
        SigmaStar = Kxx - kall / (Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) * kall';
        
                
        ypred2 = kall /(Sigma + eye(size(GammaMap))/(GammaMap + ridge*eye(size(Sigma)))) *fmap;
%         stdpred = (1 - sum((kall*iK) .* kall, 2)).^.5;
            stdpred2 = diag(SigmaStar).^.5;
         
        
        figure
        clf, plot(x, r)
        hold on, plot(xtest(1:i), ytest(1:i), 'r*')
        hold on, plot(xtest(ixRew), ytest(ixRew), 'ko')
        
     
        hold on, plot(x, ypred2, 'k');
        hold on, plot(x, ypred2 + 2*stdpred2, 'k--');
        hold on, plot(x, ypred2 - 2*stdpred2, 'k--');
        
            title('Gridsearch')
        
        
       
    end
end




% 
% ridge = 1e-10;
% profile on
% [f, df] = pref_loghyp_opt(loghyp(1:3), xtest(:), prefs, fPrior, ridge, .2);
% profile viewer
% objfun = @(lh) pref_loghyp_opt(lh, xtest(:), prefs, fPrior, ridge, .2);
% numericalGradient(objfun, loghyp(1:3))

% disp(['Minimal function values (fminunc/grid): ', num2str([fopt , min(min(min(f)))])])
% disp(['optimal solutions sigP (fminunc/grid): ', num2str([exp(logsigopt(1)) sigpOpt])])
% disp(['optimal solutions sigA (fminunc/grid): ', num2str([exp(logsigopt(2)) sigaOpt])])
% disp(['optimal quantile (fminunc/grid): ', num2str([exp(logsigopt(3)) quantileOpt])])
% disp(['Optimization times (fminumnc/grid)" ', num2str([fminuncTime, gridTime])])
