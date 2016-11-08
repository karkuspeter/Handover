function lgp = localGPUpdate(lgp, x, y, wgen, sizeLimit, onlineTrain, hypInit)
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
%	wgen: activation limit to generate new local model [1] in [0-1]
%   sizeLimit: maximum number of datapoints contained in a local model [1]
%   onlineTrain: train hyperparameters online, will be very slow, but good
%       if you do a first experiment [1]
%   hypInit: the initial log hyperparamters for a new local model
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
            knew = exp(2*lgp{j}.loghyp(d+1, imax))*exp(-.5*maha(x, lgp{j}.X{imax}, W));
            if size(lgp{j}.X{imax}, 1) > sizeLimit
                randIx = ceil(rand(1)*sizeLimit);
                lgp{j}.X{imax} = [lgp{j}.X{imax}(1:(randIx-1), :); lgp{j}.X{imax}((randIx+1):end, :)];
                lgp{j}.Y{imax} = [lgp{j}.Y{imax}(1:(randIx-1), :); lgp{j}.Y{imax}((randIx+1):end, :)];
                dummy = invupdatered(lgp{j}.iK{imax}, randIx);
                lgp{j}.iK{imax} = invupdateapp(dummy, knew([1:(randIx-1), (randIx+1):sizeLimit])',...
                    knew([1:(randIx-1), (randIx+1):sizeLimit]), knew(end)+(exp(2*lgp{j}.loghyp(d+2, imax)) + 1e-6));
            else
                lgp{j}.iK{imax} = invupdateapp(lgp{j}.iK{imax}, knew(1:end-1)', knew(1:end-1), knew(end)+(exp(2*lgp{j}.loghyp(d+2, imax)) + 1e-6));
            end
            n = size(lgp{j}.X{imax}, 1);
            
            % checking inverse accuracy
%             L = chol(exp(2*lgp{j}.loghyp(d+1, imax))*exp(-.5*maha(lgp{j}.X{imax}, lgp{j}.X{imax}, W)) + (exp(2*lgp{j}.loghyp(d+2, imax)) + 1e-6) * eye(n)); 
%             norm(lgp{j}.iK{imax} - eye(n)/L/L')
			lgp{j}.c{imax} = mean(lgp{j}.X{imax}, 1);
			% ------ Train hyperparameters here ------
            if onlineTrain
                if lgp{j}.change{imax} == 10 || n == 5
                    if n == 5
                        hypLocal = getFullGPModel(lgp{j}.X{imax}, lgp{j}.Y{imax}, round(400/n));
                    else %reuse
                        hypLocal = getFullGPModel(lgp{j}.X{imax}, lgp{j}.Y{imax}, round(400/n), lgp{j}.loghyp(:, imax));
                    end
                    lgp{j}.loghyp(:, imax) = hypLocal;
                    lgp{j}.change{imax} = 0;
                    W = diag(exp(lgp{j}.loghyp(1:d, imax)).^-2);
                    L = chol(exp(2*lgp{j}.loghyp(d+1, imax))*exp(-.5*maha(lgp{j}.X{imax}, lgp{j}.X{imax}, W)) + (exp(2*lgp{j}.loghyp(d+2, imax)) + 1e-6) * eye(n)); 
                    lgp{j}.iK{imax} = eye(n)/L/L';
                end
            end
			% ----------------------------------------
			lgp{j}.alpha{imax} = lgp{j}.iK{imax}*lgp{j}.Y{imax};
            lgp{j}.change{imax} = lgp{j}.change{imax} + 1;
		else % create new local model
			lgp{j}.c{nolgp + 1} = x;
			lgp{j}.X{nolgp + 1} = x;
			lgp{j}.Y{nolgp + 1} = y(j);
            lgp{j}.loghyp(:, nolgp + 1) = hypInit(:, j);
			lgp{j}.iK{nolgp + 1} = 1/(1e-6 + exp(2*lgp{j}.loghyp(end-1, nolgp + 1)) + exp(2*lgp{j}.loghyp(end, nolgp + 1)));
			lgp{j}.alpha{nolgp + 1} = 0;
            lgp{j}.change{nolgp + 1} = 1;
		end
	else
		lgp{j}.c{nolgp + 1} = x;
		lgp{j}.X{nolgp + 1} = x;
		lgp{j}.Y{nolgp + 1} = y(j);
		lgp{j}.loghyp(:, nolgp + 1) = hypInit(:, j);
		lgp{j}.iK{nolgp + 1} = 1/(1e-6 + exp(2*lgp{j}.loghyp(end-1, nolgp + 1)) + exp(2*lgp{j}.loghyp(end, nolgp + 1)));
		lgp{j}.alpha{nolgp + 1} = 0;
        lgp{j}.change{nolgp + 1} = 1;
	end


end