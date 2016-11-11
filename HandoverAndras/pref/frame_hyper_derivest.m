clear all, close all

rew = @(x) 2*(.5*exp(- (x - 2) * 1 * (x -2)) + exp(- x * 4 * x) -0.1174);
addpath('../gp')
addpath('..')
addpath('../DERIVESTsuite/')

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
N =50;
for i = 1:N
    
    ixRand = randperm(length(ytest), 2);
    
    if (ytest(ixRand(1)) ) >= (ytest(ixRand(2)))
        prefs(i, :) = ixRand;
    else
        prefs(i, :) = [ixRand(2) ixRand(1)];
    end
    
end
ridge = 1e-4;



for i = sampleNum
    
    
    ixRew = 1:5:i;
    fPrior = [ixRew(:), ytest(ixRew(:))'];
    
    %           fPrior = [];
    %           ixRew = [];
    
    w = medianTrick(xtest', 0.1);
    loghyp = log([1, .5, w]);
    
    options = optimoptions('fminunc', 'Algorithm','trust-region','GradObj','on','Hessian', 'off', 'MaxFunEvals', 100, 'TolX', 1e-3, 'TolFun', 1e-2);
    
    
    %     [f,df, H ] = pref_loghyp_derivest(loghyp, xtest', prefs, fPrior, ridge);
    %     [f1,df1 ] = pref_loghyp_derivest_maxStep(loghyp, xtest', prefs, fPrior, ridge, 1);
    %     [f2,df2 ] = pref_loghyp_numGrad(loghyp, xtest', prefs, fPrior, ridge, 1);
    for k = 1:2
        if k == 1
            optfun = @(lh) pref_loghyp_derivest_maxStep(lh, xtest', prefs, fPrior, ridge, 1);
        else
            optfun = @(lh) pref_loghyp_numGrad(lh, xtest', prefs, fPrior, ridge, 1);
        end
        tic
        [loghyp_opt, fopt, ~, optimOutput] = fminunc(optfun, loghyp, options);
        timeTaken = toc;
        optimOutput
        
        W = exp(loghyp_opt(3))^-2;
        sig = exp(loghyp_opt(1));
        sigma2 = exp(loghyp_opt(2));
        
        Sigma =  exp(-.5 * maha(xtest(1:i)', xtest(1:i)', W)) ;
        Sigma = Sigma + eye(size(Sigma)) * ridge;
        
        iK = eye(size(Sigma))/(Sigma);
        
        [fmap, ~, GammaMap] = nr_plgp_wPrior(zeros(length(xtest), 1), prefs, Sigma, sig, fPrior, sigma2);
        
        
        kall = exp(-.5 * maha(x(:), xtest(1:i)', W));
        Kxx =exp(-.5 * maha(x(:), x(:), W)) ;
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
        
        if k == 1
            title(['Derivest, fopt = ', num2str(fopt), ', time taken: ', num2str(timeTaken), ' sec'])
            disp('Derivest loghyp_opt')
            loghyp_opt
        else
            title(['NumGrad, fopt = ', num2str(fopt), ', time taken: ', num2str(timeTaken), ' sec'])
            disp('NumGrad loghyp_opt')
            loghyp_opt
        end
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
