function [f,df, H ] = pref_loghyp_derivest_maxStep(loghyp, x, prefs, fPrior, ridge, maxStep)

if ~isempty(fPrior)
    absoluteFeedback = 1;
else
    absoluteFeedback = 0;
end

sigp = exp(loghyp(1));
siga = exp(loghyp(2));

w = exp(loghyp(3:end));
W = diag(w.^-2);

Sigma = exp(-.5 * maha(x, x, W)) ;
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
    
    z = (fmap(prefs(:, 1)) - fmap(prefs(:, 2)) )/sqrt(2)/sigp;
        
    cdfz = normcdf(z, 0, 1);
    
    if absoluteFeedback
        targetF = (fmap(fPrior(:, 1)) - fPrior(:, 2));
    else
        targetF = 0;
    end
    S = - sum(log(cdfz)) + .5 * fmap' * iK * fmap + .5*siga^-2 * targetF(:)' * targetF(:);%+ length(targetF)*log(siga) - .5*log(det(iK)) ;
    
    f = -S  -.5 * log( det(I + K * GammaMap));
    f = -f;
catch
    keyboard
    f= NaN;
end

% disp('function value computed')

if nargout > 1
    optfun = @(lh) pref_loghyp_derivest_maxStep(lh, x, prefs, fPrior, ridge, maxStep);
%     disp('computing gradient')
    df = gradest_MaxStep(optfun, loghyp, maxStep);
%     disp('gradient computed')
end

if nargout > 2
    optfun = @(lh) pref_loghyp_derivest_maxStep(lh, x, prefs, fPrior, ridge, maxStep);
    disp('computing Hessian')
    H = hessian_maxStep(optfun, loghyp,maxStep);
    disp('Hessian computed')
end
    