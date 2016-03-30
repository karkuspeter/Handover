function [S, dS, ddS] = grad_nr_plgp_wPrior(f, prefs, K, sig, fPrior, sigma2)
% in prefs, the index of the preferred solution is the first, the dominated
% is the second.
% in prefs, the index of the preferred solution is the first, the dominated
% is the second.


iK = eye(size(K))/(K);

for i = 1:size(prefs, 1)
    
%     z(i) = (f(prefs(i, 1)) - fPrior(prefs(i, 1)) - f(prefs(i, 2)) + fPrior(prefs(i, 2)))/sqrt(2)/sig;
    z(i) = (f(prefs(i, 1)) - f(prefs(i, 2)) )/sqrt(2)/sig;
    
end

cdfz = normcdf(z, 0, 1);

targetF = (f(fPrior(:, 1)) - fPrior(:, 2));
S = - sum(log(cdfz)) + .5 * f' * iK * f + sigma2^-2/2 * targetF(:) * targetF(:)';

dS = zeros(size(K, 1), 1);
ddS = zeros(size(K));
beta = zeros(size(K, 1), 1);

dummy = iK * f;

for i = 1:length(f)
    
    ixHitAll = bsxfun(@eq, prefs, [i i]);
    ixHit = find(any(bsxfun(@eq, prefs, [i i]), 2));
    
    dCdf = 0;
    for j = 1:length(ixHit)
        preferred = 2* (ixHitAll(ixHit(j), 1) == 1) -1;
        dCdf = dCdf - preferred * normpdf(z(ixHit(j)), 0, 1)/cdfz(ixHit(j))/sqrt(2)/sig;
    end
    
    dS(i) =  dCdf + dummy(i);
    
    ixTarget = find(fPrior(:, 1) == i);
    if ~isempty(ixTarget)
        dS(i) = dS(i) + sigma2^-2 * (f(i) - fPrior(ixTarget, 2));
    end
    
    beta(i) = -dCdf;
    
    
    for j = 1:length(f)
        
        ddCdf = 0;
        ixHitAll = bsxfun(@eq, prefs, [i i]);
        ixHit = find(any(ixHitAll, 2));
        
        if i == j
            for k = 1:length(ixHit)
                ddCdf = ddCdf + .5 / sig^2 * (normpdf(z(ixHit(k)), 0, 1)^2/cdfz(ixHit(k))^2 + z(ixHit(k)) * normpdf(z(ixHit(k)), 0, 1)/cdfz(ixHit(k)));
            end
        end
        
        ixHitAll2 = or(bsxfun(@eq, prefs, [i, j]), bsxfun(@eq, prefs, [j, i]));
        ixHit2 = find(all(ixHitAll2, 2));
        for k = 1:length(ixHit2)
            loc_ddCdf = - .5 / sig^2 * (normpdf(z(ixHit2(k)), 0, 1)^2/cdfz(ixHit2(k))^2 + z(ixHit2(k)) * normpdf(z(ixHit2(k)), 0, 1)/cdfz(ixHit2(k)));
            %                 if isnan(loc_ddCdf)
            %                     keyboard
            %                 end
            ddCdf = ddCdf + loc_ddCdf;
        end
        ddS(i, j) = iK(i, j) + ddCdf;
        
        if i == j
            ixTarget = find(fPrior(:, 1) == i);
            if ~isempty(ixTarget)
                ddS(i, j) = ddS(i, j) + sigma2^-2;
            end
        end
        
        %             if isnan(ddS(i, j))
        %                 keyboard
        %             end
    end
   
end
