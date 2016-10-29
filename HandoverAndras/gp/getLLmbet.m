function [L, Lm, bet] = getLLmbet(x, y, px, hyp)
% 
% Precalculates some important parameters for 
% prediction with pseudo inputs.
%
% Inputs:
%   x: training input [n, d]
%   y: training target [n, 1]
%   px: the pseudo inputs [m, d]
%   hyp: log hyper parameters (b, c, sig) [d+2, 1]
%
% Outputs:
%   L: [m, m]
%   Lm: [m, m]
%   bet: [m, 1]

[n, d] = size(x);
[m, d] = size(px);
sigf_sq = exp(hyp(end-1));
sign_sq = exp(hyp(end));
K = kern(px, px, hyp) + 1e-6*eye(m);
L = chol(K)';
K = kern(px, x, hyp);
V = L\K;
ep = 1 + (sigf_sq*ones(n, 1)-sum(V.^2,1)')/sign_sq;
V = V./repmat(sqrt(ep)', m, 1); y = y./sqrt(ep);
Lm = chol(sign_sq*eye(m) + V*V')';
bet = Lm\(V*y);