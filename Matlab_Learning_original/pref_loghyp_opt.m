function [f, df] = pref_loghyp_opt(loghyp, x, prefs, fPrior, ridge, siga)

if ~isempty(fPrior)
    absoluteFeedback = 1;
else
    absoluteFeedback = 0;
end


hyp = exp(loghyp);
w = hyp(1:end-2);
sigf = hyp(end-1);
sigp = hyp(end);
% siga = hyp(end);

N = size(x, 1);

n = length(w);
W = diag(w.^-2);

Sigma = sigf^2 * exp(-.5 * maha(x, x, W)) ;
K = Sigma;
I = eye(size(Sigma));
Sigma = Sigma + I *N* ridge;

if absoluteFeedback
   fdummy = (rand(length(x), 1)*(max(fPrior(:, 2))-min(fPrior(:, 2)))+min(fPrior(:, 2)));
else
   fdummy = randn(length(x), 1);
end
[fmap, ~, GammaMap] = nr_plgp_wPrior(fdummy, prefs, Sigma, sigp, fPrior, siga);

iK = I/(Sigma);

for i = 1:size(prefs, 1)
    z(i) = (fmap(prefs(i, 1)) - fmap(prefs(i, 2)) )/sqrt(2)/sigp;
end

cdfz = normcdf(z, 0, 1);

if absoluteFeedback
    targetF = (fmap(fPrior(:, 1)) - fPrior(:, 2));
else
    targetF = 0;
end
S = - sum(log(cdfz)) + .5 * fmap' * iK * fmap + .5*siga^-2 * targetF(:)' * targetF(:);

f = -S  -.5 * log( det(I + K * GammaMap));

f = -f; % for minimization algorithms

if nargout > 1
    
    dummy = I/(iK + GammaMap + I * N* ridge);
    alpha = fmap'*iK;
    % kernel function gradients
    for i = 1:length(w)
        Wloc = zeros(size(W));
        Wloc(i, i) = w^-3;
        dKdw = K .* maha(x, x, Wloc);
        df(i) = w(i)/2 * alpha * dKdw * alpha' - ...
            w(i)/2 * trace(dummy * iK * dKdw * GammaMap);
    end
    dKdsigf = 2/sigf * K;
    df(n+1) = sigf/2 * alpha * dKdsigf * alpha' - ...
            sigf/2 * trace(dummy * iK * dKdsigf * GammaMap);
       
    df = -df;
    % preference noise gradient
%     dLogPhidSigP = - (normpdf(z)./normcdf(z) .* z) / sigp ;
%     dGammaMapdSigP = -2/sigp * GammaMap;
%     df(n+2) = sigp * sum(dLogPhidSigP) - sigp /2 * trace(dummy * dGammaMapdSigP);

    objfun = @(lsp) pref_loghyp_opt([loghyp(1:end-1), lsp], x, prefs, fPrior, ridge, siga);
    
    df(n+2) = numericalGradient(objfun, log(sigp));
    
  
    
end




        
