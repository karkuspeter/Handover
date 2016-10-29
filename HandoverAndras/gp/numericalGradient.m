function [df] = numericalGradient(func, init_par)
%
% Calculates the numerical gradient df/dpar
%
% Input: 
%   func: callable function that returns a real number, input is a
%   parameter vector [func]
%   init_par: df/dpar is calculated at init_par [n, 1]
%
% Output:
%   df: the gradient [n, 1]

testRet = func(init_par);
lenRet = length(testRet);

n = length(init_par);
eps_default = 1e-4;

df = zeros(n, lenRet);

init_par = init_par(:);

for i = 1:n
    eps = eps_default;
     
        perturbation = zeros(n, 2);
        perturbation(i, :) = [-eps, eps];
        f_neg = func(init_par + perturbation(:, 1));
        f_pos = func(init_par + perturbation(:, 2));
    
        df(i, :) = (f_pos(:)' - f_neg(:)')/2/eps;

end
 