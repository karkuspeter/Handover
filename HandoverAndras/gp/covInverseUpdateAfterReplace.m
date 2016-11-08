function A = covInverseUpdateAfterReplace(iKold, Kold, knew, j)
%
% Updates the inverse of a quadratic positive semidefinite matrix when
% replacing the j-th row/column.
%
% Inputs:
%   iKold: inverse of old matrix [n, n]
%   Kold:  the old matrix [n, n]
%   knew: the new row/column [n, 1]
%   j: the index of row/column to replace [1]
%
% Outputs:
%   A: the new inverse matrix [n, n]

n = size(iKold, 1);

r = knew - Kold(j, :)';
Astar = iKold - (iKold * r *iKold(j, :)) / (1 + r'* iKold(j, :)');
A = Astar - (Astar(:, j) * r' * Astar)/(1 + r'*Astar(:, j));
% A = Astar - (Astar(j, :)' * r' * Astar)/(1 + r'*Astar(j, :)');