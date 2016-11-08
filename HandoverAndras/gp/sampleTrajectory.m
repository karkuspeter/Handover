function [traj] = sampleTrajectory(s0, L, Lm, bet, px, hyp, polpar, steps);
%
% Samples a trajectory according to the GP model
%
% Inputs:
%   s0: starting state [e, 1]
%   L: parameter for SPGP prediction [m, m, e]
%   Lm: parameter for SPGP prediction [m, m, e]
%   bet: parameter for SPGP prediction [m, e]
%   px: pseudo inputs [m, d]
%   hyp: hyper parameters (b, c, sig) [d+2, e]
%   polpar: CURRENTLY LINEAR policy parameters [e, 1]
%   steps: number of steps to sample
%
% Output:
%   traj: the sampled trajectory (states, actions) [steps, d]

s = s0;
e = length(s);
udim = size(px, 2) - e;
sign_sq = exp(hyp(end, :));
sigf_sq = exp(hyp(end-1, :));
traj = zeros(steps, e+udim);
x = [s, sum(s(:).*polpar(:))];
traj(1, :) = x;
mu = zeros(1, e);
s2 = zeros(1, e);

for i = 2:steps
    for j = 1:e
        K = kern(px, x, hyp(:, j));
        lst = L(:, :, j)\K;
        clear K
        lmst = Lm(:, :, j)\lst;
        mu(j) = sum(bet(:, j).*lmst);
        s2(j) = sigf_sq(j) +sign_sq(j) - sum(lst.^2) + sign_sq(j)*sum(lmst.^2);
    end
    s = s + mu + randn(1, e).*(s2.^.5);
    x = [s, sum(s(:).*polpar(:))];
    traj(i, :) = x;
end