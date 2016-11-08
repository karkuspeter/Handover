function [logl, ll] = validateGpModel(Hyp, input, output, verb)
% Validates a GP model by calculating the log likelihood of training targets
%
% Inputs:
%	Hyp: log hyperparameters in all target dimensions [E, D+2]
%	input: the training inputs [n, D]
%	output: the training targets [n, E]
%	verb: verbosity {0, 1}
%
% Outputs:
%	logl: the loglikelihood in all dimensions [1, E]
%	ll: the likelihood value for each output [n, E]

[n, D] = size(input);
[n, E] = size(output);

logl = zeros(1, E);
ll = zeros(n, E);

for i = 1:E
	w = exp(Hyp(i, 1:D));
	sigf = exp(Hyp(i, D+1));
	signs = exp(Hyp(i, D+2));
	K = sigf^2*exp(-.5*maha(input, input, diag(w.^-2)));
	ypred = zeros(n, 1);
	yvar = zeros(n, 1);
	iK = inv(K + signs^2 * eye(n));
	dummy = iK*output(:, i);
	for j = 1:n
		k = K(:, j); 
		ypred(j) = k'*dummy;
		yvar(j) = sigf^2 - k'*iK*k;
	end
	logl(i) = sum(-.5*(ypred-output(:, i)).^2./yvar) - sum(log(yvar.^5)) - n/2*log(2*pi);
	if verb
		disp(['Model Validation dimension #',num2str(i),'\n'])
		disp('mu-y*, var, ll')
		disp([ypred-output(:, i), yvar, exp(-.5*(ypred-output(:, i)).^2./yvar)])
		disp(['loglikelihood: ',num2str(logl(i))])
	end
	ll(:, i) = (2*pi)^(-.5)*(1./yvar.^.5).*exp(-.5*(ypred-output(:, i)).^2./yvar);
end
