function [mt, covt] = getTorqueCost(m, c, Q)
%
% Calcualtes the expected cummulated torque square cost with its covariance matrix
%
% Input:
%	m: the mean of the torque trajectory [n, d]
%	c: the covariance of the torque trajectory [d, d, n]	
%	Q: quadratic weighting matrix such that cost(t) = torque'*Q*torque [d, d]
%
% Output:
%	mt: the mean cummulated torque square [d, 1]
%	covt: the variance of the cummulated torque square [d, d]

[n, d] = size(m);

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
