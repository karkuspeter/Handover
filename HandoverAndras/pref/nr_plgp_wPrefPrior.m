function [fmap, ddS, Gamma] = nr_plgp_wPrefPrior(f, prefs, K, sig, fPrior, sigma2)
% in prefs, the index of the preferred solution is the first, the dominated
% is the second.
% in prefs, the index of the preferred solution is the first, the dominated
% is the second.

maxiter = 100;
fprev = f - 1000000;
iK = eye(size(K))/(K);
step_counter = 0;

if isempty(fPrior)
    absoluteFeedback = 0;
else
    absoluteFeedback = 1;
end

while and(norm(f-fprev) > 1e-3, step_counter < maxiter)
    for i = 1:size(prefs, 1)
        
        %     z(i) = (f(prefs(i, 1)) - fPrior(prefs(i, 1)) - f(prefs(i, 2)) + fPrior(prefs(i, 2)))/sqrt(2)/sig;
        z(i) = (f(prefs(i, 1)) - f(prefs(i, 2)) )/sqrt(2)/sig;
        
    end
    
    cdfz = normcdf(z, 0, 1);
    pdfz = normpdf(z, 0, 1);
    
    
    if absoluteFeedback
        targetF = (f(fPrior(:, 1)) - fPrior(:, 2));
    else
        targetF = 0;
    end
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
            dCdf = dCdf - preferred * pdfz(ixHit(j))/cdfz(ixHit(j))/sqrt(2)/sig;
        end
        
        dS(i) =  dCdf + dummy(i);
        
        if absoluteFeedback
            ixTarget = find(fPrior(:, 1) == i);
            if ~isempty(ixTarget)
                dS(i) = dS(i) + sigma2^-2 * (f(i) - fPrior(ixTarget, 2));
            end
        end
        
        beta(i) = -dCdf;
        
        ixHitAll = bsxfun(@eq, prefs, [i i]);
        ixHit = find(any(ixHitAll, 2));
        for j = 1:length(f)
            
            if 0
                ddS(i, j) = ddS(j, j);
                Gamma(i, j) = Gamma(j, i);
            else
                ddCdf = 0;
                
                
                if i == j
                    for k = 1:length(ixHit)
                        ddCdf = ddCdf + .5 / sig^2 * (pdfz(ixHit(k))^2/cdfz(ixHit(k))^2 + z(ixHit(k)) * pdfz(ixHit(k))/cdfz(ixHit(k)));
                    end
                end
                
                if j >=i
                    ixHitAll2 = or(bsxfun(@eq, prefs, [i, j]), bsxfun(@eq, prefs, [j, i]));
                    ixHit2 = find(all(ixHitAll2, 2));
                    ixHit2Save{i, j} = ixHit2;
                else
                    ixHit2 = ixHit2Save{j, i};
                end
                
                for k = 1:length(ixHit2)
                    loc_ddCdf = - .5 / sig^2 * (pdfz(ixHit2(k))^2/cdfz(ixHit2(k))^2 + z(ixHit2(k)) * pdfz(ixHit2(k))/cdfz(ixHit2(k)));
                    %                 if isnan(loc_ddCdf)
                    %                     keyboard
                    %                 end
                    ddCdf = ddCdf + loc_ddCdf;
                end
                if absoluteFeedback
                    if i == j
                        ixTarget = find(fPrior(:, 1) == i);
                        if ~isempty(ixTarget)
                            ddCdf = ddCdf + sigma2^-2;
                        end
                    end
                end
                Gamma(i, j) = ddCdf;
                
                if isnan(ddCdf)
                    warning('NaN in ddCdf')
                end
                
                ddS(i, j) = iK(i, j) + ddCdf;
                
            end
            
        end
        
    end
    
    %     nhess = num_hess(optfun, f, 1e-6);
    %     ngrad = num_grad(optfun, f, 1e-4);
    
    fanal = f - ddS\dS;
    %     fnum = f - nhess\dS;
    fprev = f;
    f = fanal; %fnum;
    
    
    step_counter = step_counter + 1;
    
    %     figure, plot([dS, ngrad]);
    %     figure,plotMatrix((ddS - nhess).^2);
    %
    %     figure, plot([fanal, fnum])
    %     drawnow
    %     keyboard
    
end
fmap = f;

if step_counter > maxiter
    warning('f optimization run over')
    
end

