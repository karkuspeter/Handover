function y = localGPVarSparsePredict(lgp, x, m)
%
% Predicts output given local GP model. For details see 'localGPUpdate.m'!
%
% Inputs:
%   lgp: local gp update, see 'localGPUpdate.m'
%     x: query input [1, d]
%     m: maximum number of models used for prediction [1]
% Output:
%   y: weighted prediction [1, e];

e = length(lgp);
d = length(x);
y = zeros(1, e);

for j = 1:e
    
    try
        nolgp = length(lgp{j}.X); % number of local models
    catch
        nolgp = 0;
    end
    
    mloc = m;
    mloc = min(mloc, nolgp);
           
	if nolgp > 0 % if we have local models for the current output
        activation = zeros(nolgp, 1);
        ylocal = zeros(nolgp, 1);
		for i = 1:nolgp
			W = diag(exp(lgp{j}.loghyp(1:d, i)).^-2);
			activation(i) = exp(-.5*(x - lgp{j}.c{i})*W*(x - lgp{j}.c{i})');
            k = exp(2*lgp{j}.loghyp(d+1, i)) * exp(-.5*maha(x, lgp{j}.px{i}, W));
            ylocal(i) = k*lgp{j}.alpha{i};
        end
        
        [descActivation, ix] = sort(activation, 'descend');
        ylocal = ylocal(ix);
        y(j) = sum(descActivation(1:mloc) .* ylocal(1:mloc))/sum(descActivation(1:mloc));
    else
        y(j) = 0;
    end
end