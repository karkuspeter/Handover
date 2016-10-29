function [f] = pref_loghyp_optCV(loghyp, x, prefs, fPrior, ridge, siga, samples)

if ~isempty(fPrior)
    absoluteFeedback = 1;
else
    absoluteFeedback = 0;
end

CV = 5;

hyp = exp(loghyp);
w = hyp(1:end-3);
sigf = hyp(end-2);
sigp = hyp(end-1);
siga = hyp(end);

n = length(w);
W = diag(w.^-2);

Sigma = sigf^2 * exp(-.5 * maha(x, x, W)) ;
K = Sigma;
I = eye(size(Sigma));
Sigma = Sigma + I * ridge;

if absoluteFeedback
   fdummy = (rand(length(x), 1)*(max(fPrior(:, 2))-min(fPrior(:, 2)))+min(fPrior(:, 2)));
else
   fdummy = randn(length(x), 1);
end
[fmap, ~, GammaMap] = nr_plgp_wPrior(fdummy, prefs, Sigma, sigp, fPrior, siga);
N = length(fmap);
f = 0;
for i = 1:CV
    
    ixTest = i:CV:N;
    ixTrain = setdiff(1:N, ixTest);
    
    kall = sigf^2 *exp(-.5 * maha(x(ixTest, :), x(ixTrain, :), W))';
        Kxx = sigf^2 *exp(-.5 * maha(x(ixTest, :), x(ixTest, :), W)) ;
        SigmaStar = Kxx - kall' / (Sigma(ixTrain, ixTrain) + eye(size(GammaMap(ixTrain, ixTrain)))/(GammaMap(ixTrain, ixTrain) + ridge*eye(size(Sigma(ixTrain, ixTrain))))) * kall;
        ypred = (kall' / Sigma(ixTrain, ixTrain)) * (fmap(ixTrain) );
        
        SigmaStar = (SigmaStar + SigmaStar')/2;
        
        mu0 = ypred;
        Sigma0 = SigmaStar;
        
        mu1 = fmap(ixTest);
        Sigma1 = (Sigma(ixTest, ixTest) + eye(size(GammaMap(ixTest, ixTest)))/(GammaMap(ixTest, ixTest) + ridge*eye(size(Sigma(ixTest, ixTest)))));

        % kl of mu0,Sigma0 || mu1, Sigma1
        
        kl = .5 * (trace(Sigma1\Sigma0) + (mu1-mu0)'/Sigma1*(mu1-mu0) - length(w) + log(det(Sigma1)/det(Sigma0)));
        ikl = .5 * (trace(Sigma0\Sigma1) + (mu1-mu0)'/Sigma0*(mu1-mu0) - length(w) + log(det(Sigma0)/det(Sigma1)));
        
        kl = (mu1-mu0)'*(mu1-mu0);
        f = f+kl;
        
end
        
