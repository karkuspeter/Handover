% Demonstrating how the exponential scaling changes the weight distribution
% based on the reward distribution (uniform, gaussian, etc.)
%
% w = exp(r/eta);
%

clear all, close all

% uniform reward distribution
r = 0:.1:10; rlen = length(r);
eta = [1:.5:10]'; etalen = length(eta);

w = exp(bsxfun(@rdivide, r, eta));

rr = repmat(r, etalen, 1);
eeta = repmat(eta, 1, rlen);

surf(rr, eeta, log(w))

xlabel('reward'), ylabel('\eta'), zlabel('log w')
title('Uniform reward distro')

% uniform reward distribution
x = 0:.1:10; r = 5*exp(-.5/1.5^2 * (x - 5).^2);
rlen = length(r);
eta = [1:.5:10]'; etalen = length(eta);

w = exp(bsxfun(@rdivide, r, eta));

rr = repmat(r, etalen, 1);
eeta = repmat(eta, 1, rlen);

figure
surf(rr, eeta, log(w))

xlabel('reward'), ylabel('\eta'), zlabel('log w')




