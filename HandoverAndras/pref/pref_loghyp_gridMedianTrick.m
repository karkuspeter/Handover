function [f, sigpOpt, sigaOpt, quantileOpt] = pref_loghyp_gridMedianTrick(x, prefs, fPrior, ridge, resolutionSigP, resolutionSigA, resolutionQuantile, plote)

if ~isempty(fPrior)
    absoluteFeedback = 1;
else
    absoluteFeedback = 0;
end

sigf = 1;
if absoluteFeedback
    fdummy = (rand(length(x), 1)*(max(fPrior(:, 2))-min(fPrior(:, 2)))+min(fPrior(:, 2)));
else
    fdummy = randn(length(x), 1);
end


minOpt = Inf;

for ii = 1:length(resolutionSigP)
    sigp = resolutionSigP(ii);
    for jj = 1:length(resolutionSigA)
        siga = resolutionSigA(jj);
        for kk = 1: length(resolutionQuantile)
            w = medianTrick(x, resolutionQuantile(kk));
            W = diag(w.^-2);
            
            Sigma = sigf^2 * exp(-.5 * maha(x, x, W)) ;
            K = Sigma;
            I = eye(size(Sigma));
            Sigma = Sigma + I * ridge;
            iK = I/(Sigma);
            
            
            
            try
                [fmap, ~, GammaMap] = nr_plgp_wPrior(fdummy, prefs, Sigma, sigp, fPrior, siga);
                
                %             for i = 1:size(prefs, 1)
                z = (fmap(prefs(:, 1)) - fmap(prefs(:, 2)) )/sqrt(2)/sigp;
                %             end
                
                cdfz = normcdf(z, 0, 1);
                
                if absoluteFeedback
                    targetF = (fmap(fPrior(:, 1)) - fPrior(:, 2));
                else
                    targetF = 0;
                end
                
                S = - sum(log(cdfz)) + .5 * fmap' * iK * fmap -.5*siga^-2 * targetF(:)' * targetF(:);% - length(targetF)*log(siga)  - .5*log(det(iK)) ;
                
                f(ii, jj, kk) = -S  -.5 * log( det(I + K * GammaMap));
                
                f(ii, jj, kk) = -f(ii, jj, kk); % we are minimizing
                
                if f(ii, jj, kk) < minOpt
                    minOpt = f(ii, jj, kk);
                    sigpOpt = sigp;
                    sigaOpt = siga;
                    quantileOpt = resolutionQuantile(kk);
                end
            catch
                f(ii, jj, kk) = NaN;
            end
          
        end
    end
end
    
if plote
figure, contourf(f(:, :, find(quantileOpt == resolutionQuantile))), colorbar, xlabel('sigA'), ylabel('sigP')
h = get(gca);
for i = 1:length(resolutionSigA)
    resolutionSigmaACell{i} = num2str(resolutionSigA(i), 4);
    resolutionSigmaPCell{i} = num2str(resolutionSigP(i), 4);
end
set(gca, 'XTickLabel',  [resolutionSigmaACell']);
set(gca, 'YTickLabel',  [resolutionSigmaPCell']);
    
end
