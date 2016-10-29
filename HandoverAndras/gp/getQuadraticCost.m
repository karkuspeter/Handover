function [mt, covt] = getQuadraticCost(m, c, Q)
%
% Calcualtes the expected cummulated square cost with its covariance matrix for a trajectory.
%
% Input:
%	m: the mean of the trajectory [n, d]
%	c: the covariance of the trajectory [d, d, n]	
%	Q: quadratic weighting matrix such that cost(t) = x'*Q*x, [d, d]
%
% Output:
%	mt: the mean cummulated square [1]
%	covt: the variance of the cummulated square [1]

[n, d] = size(m);

if d == 1 % in this case size(Q) = [1, n], we extend to [1, 1, n]
	dummy = ones(1, 1, n);
	dummy(1, 1, :) = c;
	c = dummy;
end

if nargin == 2
	Q = eye(d);
end

mt = 0;
covt = 0;

for i = 1:n
	mt = mt + trace(c(:, :, i)*Q) + m(i, :)*Q*m(i, :)';
	vart = trace(2*Q*c(:, :, i)*Q*c(:, :, i)) + 4*m(i, :)*Q*c(:, :, i)*Q*m(i, :)';
	covt = covt + vart;
end
