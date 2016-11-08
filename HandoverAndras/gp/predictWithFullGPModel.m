function [M, V, K] = predictWithFullGPModel(loghyp, x, y, z)
%
% Predict the mean and standard deviation at test input z using the
% optimized hyperparameters and the training data.
%
% Input:
%   Hyp: a struct of log hyperaparameters, see getFullGPModel.m for details
%   x: training input [n, d]
%   y: training target [n, e]
%   z: test input [m, d]


[n, e] = size(y);
hyp = exp(loghyp);

for i = 1:e
    w = hyp(1:end-2, i);
    W = diag(w.^-2);
    sigf = hyp(end-1, i);
    sign = hyp(end, i);
    
    K = sigf^2 * exp(-.5 * maha(x, x, W));
    Ky = K + sign^2*eye(size(K));
    k = sigf^2 * exp(-.5 * maha(x, z, W));
    
    M(:, i) = k'/Ky * y(:, i);
    V(:, i) = sigf^2 - sum(bsxfun(@times, k'/Ky, k'), 2);
    
end

% 
% d = size(x, 2);
% e = size(y, 2);
% m = size(z, 1);
% covfunc = @covSEard; 
% likfunc = @likGauss; 
% n = size(x, 1);
% 
% M = zeros(m, e);
% S = zeros(m, e);
% xx = [x; z];
% 
% for i = 1:e
%     hyp.cov = hypIn(1:end-1)';
%     hyp.lik = hypIn(end);
%     hyp.mean = [];
%     Hyp{i} = hyp;
%     
%     [m, s2] = gp(Hyp{i}, @infExact, [], covfunc, likfunc, x, y(:, i), z);
%     M(:, i) = m;
%     S(:, i) = sqrt(s2);
%     
%     b = (Hyp{i}.cov(1:end-1)).^-2;
%     c = (Hyp{i}.cov(end)).^2;
%     
% 	if nargout == 3
% 		for j = 1:size(xx,1)
% 			for k = j:size(xx, 1)
% 				dummy = (xx(j, :)-xx(k, :));
% 				K(j, k, i) = exp(-.5*dummy*diag(b)*dummy');
% 				K(k, j, i) = K(j, k, i);
% 			end
% 		end
%     end
% end
