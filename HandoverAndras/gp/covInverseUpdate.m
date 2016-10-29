function invNew = covInverseUpdate(invOld, colNew, diagNew);
%
% Efficient matrix inverse update for symmetric matrices, when adding a new element to the last col/row, such that invNew = [invOld^-1, colNew; colNew^T, diagNew]^-1
%
% Inputs:
%	invOld: the old inverse Matrix [n, n]
%	colNew: the new column in the [n, 1]
%	diagNew: the new variance [1, 1]
%
% Output:
%   invNew: the new inverse covariance [n+1, n+1]

alpha = invOld*colNew;
gamma = diagNew - colNew'*alpha;
invNew = 1/gamma* [gamma*invOld + alpha*alpha', -alpha
					-alpha', 1];

