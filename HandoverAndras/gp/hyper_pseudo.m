function [hyp, hyp_orig, px] = hyper_pseudo(x, y, m, steps);
% Calculates the hyperparameters and pesudo inputs of a GP model
% (Requires Ed Snelson's 'SPGP_dist' library, assumed to be in the path)
% 
% Input:
%	x: training samples [n, d]
%   y: target samples [n, e]
%	m: number of pseudo samples
%
% Output:
%	hyp: log hyperparameters [d+2, e]
%   hyp_orig: log hyperparameters (b, c, sig) [d+2, e]
%   px: the pseudo inputs [m, d, e]

[n, d] = size(x);
[n, e] = size(y);

% initialize pseudo-inputs to a random subset of training inputs
randix = randperm(n);
xb_init = x(randix(1:m), :);

% initialize hyperparameters sensibly (see spgp_lik for how
% the hyperparameters are encoded)
hyp = zeros(d+2, e);
hyp_orig = zeros(d+2, e);
px = zeros(m, d, e);
for i = 1:e
	hyp_init = zeros(d+2, 1);
	hyp_init(1:d) = -2*log(std(x)'); % log 1/(lengthscales)^2
	hyp_init(d+1) = log(std(y(:, i))); % log size 
	hyp_init(d+2) = log(std(y(:, i))/10); % log noise

	% optimize hyperparameters
	w_init = [reshape(xb_init, m*d, 1); hyp_init];
	[w, f] = minimize(w_init, 'spgp_lik', -steps, y(:, i), x, m);
	xb = reshape(w(1:m*d, 1), m, d);
	px(:, :, i) = xb;
    hyp_orig(:, i) = w(m*d+1:end, 1);
	hyp(:, i) = w(m*d+1:end, 1);
	hyp(1:d, i) = log(exp(hyp(1:d, i)).^-.5); % b --> w
	hyp(d+1, i) = log(exp(hyp(d+1, i)).^.5); % c --> sig_f
	hyp(d+2, i) = log(exp(hyp(d+2, i)).^.5); % sig --> sig_n
end

