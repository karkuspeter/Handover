function [f ] = pref_loghyp_MedianTrick(logsigPsigA, x, prefs, fPrior, ridge)

if ~isempty(fPrior)
    absoluteFeedback = 1;
else
    absoluteFeedback = 0;
end

sigp = exp(logsigPsigA(1));
siga = exp(logsigPsigA(2));

w = medianTrick(x, exp(logsigPsigA(3)));
sigf = 1;
N = size(x, 1);
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
iK = I/(Sigma);




try
    [fmap, ~, GammaMap] = nr_plgp_wPrior(fdummy, prefs, Sigma, sigp, fPrior, siga);
    
    for i = 1:size(prefs, 1)
        z(i) = (fmap(prefs(i, 1)) - fmap(prefs(i, 2)) )/sqrt(2)/sigp;
    end
    
    cdfz = normcdf(z, 0, 1);
    
    if absoluteFeedback
        targetF = (fmap(fPrior(:, 1)) - fPrior(:, 2));
    else
        targetF = 0;
    end
    S = - sum(log(cdfz)) + .5 * fmap' * iK * fmap - .5*siga^-2 * targetF(:)' * targetF(:) + length(targetF)*log(siga)  - .5*log(det(iK)) ;
    
    f = -S  -.5 * log( det(I + K * GammaMap));
    f = -f;
catch
    f= NaN;
end





