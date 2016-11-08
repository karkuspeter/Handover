function lgp = localGPVarSparseUpdate(lgp, x, y, wgen, sizeLimit, m, steps)
%
% Updates the local GP model(s)
%
% Inputs:
%	lgp: the structure containing elements:
%		c: center [1, d]
%		loghyp: log hyperparameters [d+2, e] (w, sigf, sign)
%		iK: inverse (K+sig^2I) [n, n]
%		alpha: mean factor [n, 1]
%		X: training samples	[n, d]
%		Y: training targets [n, e]
%	x: new input [1, d]
%	y: new target [1, e]
%	wgen: activation limit to generate new local model [1] (in [0, 1])
%   sizeLimit: maximum number of training samples stored for a model [1]
%   m: number of pseudo inputs for a local model [1]
%   steps: hyperparameter optimization steps [1]
%
% Outputs:
%	lgp: updated local model

d = length(x);
e = length(y);

for j = 1:e % for each output dimension

    try
        nolgp = length(lgp{j}.X); % number of local models
    catch
        nolgp = 0;
    end
           
	if nolgp > 0 % if we have local models for the current output
		activation = zeros(nolgp, 1);
		for i = 1:nolgp
			W = diag(exp(lgp{j}.loghyp(1:d, i)).^-2);
			activation(i) = exp(-.5*(x - lgp{j}.c{i})*W*(x - lgp{j}.c{i})');
 		end
		[v, imax] = max(activation); % imax: index of maximally activated local model for the current output dimension

		if v > wgen
			lgp{j}.X{imax} = [lgp{j}.X{imax}; x];
    		lgp{j}.Y{imax} = [lgp{j}.Y{imax}; y(j)];
            if size(lgp{j}.X{imax}, 1) > sizeLimit
                lgp{j}.X{imax} = lgp{j}.X{imax}(2:end, :);
                lgp{j}.Y{imax} = lgp{j}.Y{imax}(2:end);
            end
			n = size(lgp{j}.X{imax}, 1);
			% ------ Train hyperparameters here ------
            if n == 50 || lgp{j}.change{imax} == 50 % First training
                mloc = min(n, m);
                keyboard
                [hyp, px] = hyper_variationalPseudo_commonInputs(...
                    lgp{j}.X{imax}, lgp{j}.Y{imax}, mloc, steps);
                
                lgp{j}.loghyp(:, imax) = ed2carl(hyp);
                lgp{j}.px{imax} = px;
                mf = variationalSparsePrecomputations(hyp, {px}, ...
                    lgp{j}.X{imax}, lgp{j}.Y{imax}, 0);
                lgp{j}.alpha{imax} = mf{1};
                lgp{j}.change{imax} = 0;
                
%                 L = chol(exp(2*lgp{j}.loghyp(d+1, imax))*exp(-.5*maha(lgp{j}.px{imax}, lgp{j}.px{imax}, W)) + (exp(2*lgp{j}.loghyp(d+2, imax)) + 1e-6) * eye(mloc)); 
            end
            % Update hyperparameters if the new sample has low likelihood
            lgp{j}.c{imax} = mean(lgp{j}.px{imax}, 1);
            
			% ----------------------------------------
            lgp{j}.change{imax} = lgp{j}.change{imax} + 1;
		else % create new local model
			lgp{j}.c{nolgp + 1} = x;
			lgp{j}.X{nolgp + 1} = x;
            lgp{j}.px{nolgp + 1} = x;
			lgp{j}.Y{nolgp + 1} = y(j);
			lgp{j}.loghyp(:, nolgp + 1) = log([ones(d, 1)*1e6; 1e-3; 1e-4]);
			lgp{j}.alpha{nolgp + 1} = 0;
            lgp{j}.change{nolgp + 1} = 1;
		end
	else
		lgp{j}.c{nolgp + 1} = x;
		lgp{j}.X{nolgp + 1} = x;
        lgp{j}.px{nolgp + 1} = x;
		lgp{j}.Y{nolgp + 1} = y(j);
		lgp{j}.loghyp(:, nolgp + 1) = log([ones(d, 1)*1e6; 1e-3; 1e-4]);
		lgp{j}.alpha{nolgp + 1} = 0;
        lgp{j}.change{nolgp + 1} = 1;
	end


end