function [hyp, px, hyper_traj] = hyper_pseudo_commonInputs(x, y, m, steps, hyp, px)
% Calculates the hyperparameters and pesudo inputs of a GP model
% (Requires Ed Snelson's 'SPGP_dist' library, assumed to be in the path)
% It uses the same pseudo inputs for all output dimensions.
% 
% Input:
%	x: training samples [n, d]
%   y: target samples [n, e]
%	m: number of pseudo samples
%   steps: number of maximum function evaluations
%   Additionally:
%       hyp: the initial hyperparameters [d+2, e]
%       px: the initial pseudo inputs [m, d]
%
% Output:
%	hyp: log hyperparameters [d+2, e]
%   px: the pseudo inputs [m, d]
%   hyper_traj: the trajectory of hyperparameters + pseudo inputs 
%               [steps, m*d + e*(d+2) + 1]
%                       pseudo + hyper + sum neg_log_llhood

[n, d] = size(x);
[n, e] = size(y);

if nargin < 5
    % initialize pseudo-inputs to a random subset of training inputs
    randix = randperm(n);
    xb_init = x(randix(1:m), :);
    
    % initialize hyperparameters sensibly (see spgp_lik for how
    % the hyperparameters are encoded)
    hyp = zeros(d+2, e);
    px = zeros(m, d);
    hyp_init = [];
    
    for i = 1:e
        hyp_init = [hyp_init; -2*log(std(x)'); % log 1/(lengthscales)^2
            log(std(y(:, i))^2); % log sf^2
            log((std(y(:, i))^2)/5)]; % log sn^2
    end
else
    hyp_init = reshape(hyp, e*(d+2), 1);
    xb_init = px;
end

% optimize hyperparameters
w_init = [reshape(xb_init, m*d, 1); hyp_init];

[w, f, dummy, hyper_traj] = minimize_with_progress_output(w_init, 'spgp_lik_commonPseudoInputs', -steps, y, x, m);
xb = reshape(w(1:m*d, 1), m, d);
px = xb;


for i = 1:e
	hyp(:, i) = w((m*d+1 + (i-1)*(d+2)) : (m*d + i*(d+2)));
	hyp(1:d, i) = log(exp(hyp(1:d, i)).^-.5); % b --> w
	hyp(d+1, i) = log(exp(hyp(d+1, i)).^.5); % c --> sig_f
	hyp(d+2, i) = log(exp(hyp(d+2, i)).^.5); % sig --> sig_n
end

