function [sqrew, vrew] = getQuadraticCostIndividual(m, c);
%
% Calcualtes the expected square cost with its covariance matrix for an n dimensional input.
%
% Input:
%	m: the mean [n, d]
%	c: the covariance [d, d, n]	or [d, n]. if [d, n] --> diagonal [d, d, n]
%	Q: quadratic weighting matrix such that cost(t) = x'*Q*x, [d, d]
%
% Output:
%	mt: the mean cummulated square [n, 1]
%	covt: the variance of the cummulated square [n, 1]

[n, d] = size(m);

if size(c, 1) ~= size(c, 2) % in this case size(Q) = [d, n], we extend to diagonal [d, d, n]
	dummy = ones(d, d, n);
	for i = 1:n
		dummy(:, :, i) = diag(c(:, i));
	end
	c = dummy;
end

if nargin == 2
	Q = eye(d);
end

sqrew = zeros(n, 1);
vrew = zeros(n, 1);

for i = 1:n
	sqrew(i) = trace(c(:, :, i)*Q) + m(i, :)*Q*m(i, :)';
	vrew(i) = trace(2*Q*c(:, :, i)*Q*c(:, :, i)) + 4*m(i, :)*Q*c(:, :, i)*Q*m(i, :)';
end
