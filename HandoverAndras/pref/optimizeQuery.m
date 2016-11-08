function [f, df] = optimizeQuery(xquery, x, fmap, GammaMap, loghyp, ridge)

hyp = exp(loghyp);

w = hyp(1:end-2); W = diag(w^-2);
sigf = hyp(end-1);
sign = hyp(end); 


Sigma = sigf^2 * exp(-.5 * maha(x, x, W));
I = eye(size(Sigma));
iK = I/(Sigma + I * ridge);
k = sigf^2 * exp(-.5 * maha(xquery, x, W));


MuStar = k * iK * fmap;
alpha =   (Sigma + I * ridge + I/(GammaMap + I* ridge)) \ k';
SigmaStarSq = sigf^2 - k * alpha;

zStar = (MuStar - fmap) / sqrt(2 *sign^2 + SigmaStarSq);
P = normcdf(zStar);

f = -sum(log(P))/length(P);


if nargout > 1
    
    dkstardX = -k(:) .* ( bsxfun(@minus, xquery(:)', x) * W);
    
    dMudX = dkstardX' * iK *fmap;
    dSigmaSqdX = -2 * dkstardX' * alpha;
    dummy = sqrt(2*sign^2 + SigmaStarSq);
    dDummydX = 1/dummy * dSigmaSqdX;
    
    
%     dErfdZ = 2/sqrt(pi) .* exp(-zStar.^2);
    dPdZ = 1/sqrt(2*pi) .* exp(-zStar.^2 /2);
    
    dZdX = 1/dummy * (dMudX * dummy - dDummydX * (MuStar - fmap));
    
    dPdXstar = repmat(dPdZ(:)', length(w), 1) .* dZdX';
    df = -1/length(P) * dPdXstar * (1./P(:));
    
    
    
end