function [fw, dfw] = spgp_lik_commonPseudoInputs(wext,y,x,m)
% A wrapper around 'sgpg_lik' such that for all output dimensions the
% same pseudo inputs are used. The hyperparameters will be optimized
% separately.
%
% Input: 
% 	wext: the extended hyperparameters [m*d + e*(d+2), 1]
%	y: the training targets [n, e]
%   x: tre training inputs [n, d]
%	m: the number of pseudo targets [1]
%
% Outputs:
% 	fw: the sum of negative log-likelihoods [1]
%	dfw: the gradient of wext w.r.t. its parameters [m*d + e*(d+2), 1]

[n, e] = size(y);
[n, d] = size(x);

fws = zeros(e, 1);
dfws = zeros(m*d + d+2, e);

dfw = zeros(m*d + e*(d+2), 1);
parfor i = 1:e
	wloc = [wext(1:m*d); wext((m*d+ (i-1)*(d+2)+1):(m*d+ i*(d+2)))];
	[dummy1, dummy2] = spgp_lik(wloc, y(:, i), x, m);
	fws(i) = dummy1;
	dfws(:, i) = dummy2;
end

for i = 1:e
    dfw((m*d + (i-1)*(d+2)+1):(m*d+ i*(d+2))) = dfws(m*d+1:end, i);
end

fw= sum(fws);
dfw(1:m*d) = sum(dfws(1:m*d, :), 2);

