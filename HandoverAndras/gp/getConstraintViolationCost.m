function [mt, covt] = getConstraintViolationCost(m, c, minx, maxx, Q)
% Calcualtes the expected cummulated square cost with its covariance matrix, 
% when certain constraints on the mean are violated.
%
% Input:
%	m: the mean of the trajectory [n, d]
%	c: the covariance of the trajectory [d, d, n]	
%	Q: quadratic weighting matrix such that cost(t) = x'*Q*x [d, d]
%	minx: the minimum threshold [1, d]
%	maxx: the maximum threshold [1, d]
%
% Output:
%	mt: the mean cummulated torque square [d, 1]
%	covt: the variance of the cummulated torque square [d, d]

[n, d] = size(m);

if nargin == 4
	Q = eye(d);
end
qorig = diag(Q);

mt = 0;
covt = 0;

for i = 1:n
	ix_ok = find((m(i, :) < maxx) .* (m(i, :) > minx));
	q = qorig; q(ix_ok) = zeros(1, length(ix_ok)); Q = diag(q);
	d = max([minx - m(i, :); m(i, :) - maxx]); 
	mt = mt + trace(c(:, :, i)*Q) + d*Q*d';
	vart = trace(2*Q*c(:, :, i)*Q*c(:, :, i)) + 4*d*Q*c(:, :, i)*Q*d';
	covt = covt + vart;
end
