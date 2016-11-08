function [f, df] = grad_logPref(x, loghyp, X, fmap, iKG, ridge)

sig = exp(loghyp(1));
w = exp(loghyp(2:end));
W = diag(w.^-2);

K = exp(-.5 * maha(X',X', W)) ;
K = K + eye(size(K)) * ridge;

k = exp(-.5 * maha(x(:)', X', W));
iK = eye(size(K))/(K);

fmap_star = k * iK * fmap;
var_star = 1 - k * iKG *k' + 2*sig^2;

z = bsxfun(@minus, fmap_star, fmap)./var_star.^.5;
Phi = normcdf(z, 0, 1);

f = -sum(Phi);


if nargout > 1
    
    dk = k'.* (- bsxfun(@minus, x, X)'*W);
    dvar = -2 * dk' * iKG * k(:);
    dmu = dk' * iK * fmap ;
    
    dz = 1/var_star.^.5 * dmu - (fmap_star - fmap)/2/var_star^1.5 * dvar;
    dPhi = sum(normpdf(z, 0, 1) .* dz);
    df = -dPhi;
    
end





