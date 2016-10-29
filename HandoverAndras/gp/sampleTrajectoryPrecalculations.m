function [mean_factor, var_factor] = sampleTrajectoryPrecalculations(x, y, px, hyp)
%
% Precalculations for sampling trajactory with SPGP model on GPU or CPU
%
% Inputs:
%	x: training inputs [N, d]
%	y: training trgets [N, e]
%	px: pseudo inputs  [m, d]
%	hyp: log hyper parameters (b, c, sig) [d+2, e]
%
% Outputs:
%	mean_factor: parameter for mean prediction [dim, m]	
%	var_factor: parameter for variance prediction [m, m, dim]
[N, d] = size(x);
m = size(px, 1);
dim = size(y, 2);

mean_factor = zeros(dim, m);
var_factor = zeros(m, m, dim);
for i = 1:dim
    sigf_sq = exp(hyp(end-1, i));
    sign_sq = exp(hyp(end, i));
    Km = kern(px, px, hyp(:, i)) + 1e-6*eye(m);
    iKm = inv(Km);
    Kmn = kern(px, x, hyp(:, i));
    lamb = sigf_sq*ones(N, 1) - bsxfun(@dot, Kmn, Km\Kmn)';
    Qm = Km + Kmn*(bsxfun(@times, (1./(lamb + sign_sq))', Kmn))';
    mean_factor(i, :) = (Qm\(Kmn*(y(:, i)./(lamb + sign_sq))))';
    var_factor(:, :, i) = iKm - inv(Qm+1e-6*eye(m));
end
