function Hypret = getFullGPModel(x, y, steps, hypInit)
%
% Gets the hyperparameters for non-sparsified GP model for each output
% dimensions (w, sigf, sign). For details see the GPML library.
%
% Input:
%   x: training input [n, d]
%   y: training output [n, e]
%   steps: hyperparameter optimization steps (optional) [1]
%   hyp: init of hyperparameters (optional) [d+2, e]
%
% Output:
%   Hyp: the log-hyperparameters [struct of length e]
if nargin < 3
    steps = 500;
end

[n, d] = size(x);
[n, e] = size(y);

covfunc = @covSEard; 
likfunc = @likGauss; 

Hyp = {};

if nargin < 4
    for i = 1:e
        hyp.cov = log([std(x), std(y(:, i))]);
        hyp.lik = log(std(y(:, i))/10);
        hyp.mean = [];
        Hyp{i} = hyp;
    end
else
    for i = 1:e
        hyp.cov = hypInit(1:d+1, e)';
        hyp.lik = hypInit(end, e);
        hyp.mean = [];
        Hyp{i} = hyp;
    end
end

Hypret = zeros(d+2, e);
for i = 1:e
    Hyp{i} = minimize(Hyp{i}, @gp, -steps, @infExact, [], covfunc, likfunc, x, y(:, i));
    Hypret(:, i) = [Hyp{i}.cov'; Hyp{i}.lik];
end

%     [m s2] = gp(hyp, @infExact, [], covfunc, likfunc, x', y', x');