function [dkl, idkl] = KLdivTraj(m1, c1, m2, c2)
% Calculates the KL divergence between two trajectories.
% Each state is represented with a mean and covariance.
% 
% Input:
% 	m1: mean of state [n, d]
%	c1: covariance matrix of state [d, d, n]
% 	m2: as m1
%	c2: as c1
%
% Output:
%	dkl: the KL divergence along the trajectory [n, 1]
%	idkl: the inverse KL divergence along the trajectory [n, 1]

[n, d] = size(m1);

dkl = zeros(n, 1);
idkl = zeros(n, 1);

for i = 1:n
	cov1 = c1(:, :, i);
	cov2 = c2(:, :, i);
	mn1 = m1(i, :)';
	mn2 = m2(i, :)';
	icov1 = inv(cov1);
	icov2 = inv(cov2);
	dkl(i) = .5 * (trace(icov2*cov1) + (mn2-mn1)'*icov2*(mn2-mn1) - ...
			       log(det(cov1)/det(cov2)) - d);
	idkl(i) = .5 * (trace(icov1*cov2) + (mn1-mn2)'*icov1*(mn1-mn2) - ...
			       log(det(cov2)/det(cov1)) - d);
end

