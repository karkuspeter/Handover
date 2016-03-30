function [fmap, ddS, Gamma] = nr_plgp(f, prefs, K, sig)
% in prefs, the index of the preferred solution is the first, the dominated
% is the second.
% in prefs, the index of the preferred solution is the first, the dominated
% is the second.

maxiter = 100;
fprev = f - 1000000;
iK = eye(size(K))/(K);
step_counter = 0;
while or(norm(f-fprev) > 1e-5, step_counter > maxiter)
    
    for i = 1:size(prefs, 1)
        
        z(i) = (f(prefs(i, 1)) - f(prefs(i, 2)))/sqrt(2)/sig;
%         cdfz(i) = normcdf(z(i), 0, 1, []);
        
    end
    
    cdfz = normcdf(z, 0, 1);
    
   
    S = - sum(log(cdfz)) + .5 * f' * iK * f;
    
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
            Gamma(i, j) = ddCdf;
            ddS(i, j) = iK(i, j) + ddCdf;
            
%             if isnan(ddS(i, j))
%                 keyboard
%             end
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
    disp('f optimization run over')
    keyboard
end

