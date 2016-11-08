function [xsave] = sampleTrajectoryGPU4(train_x, train_y, px, hyp, ...
    start_states, polpar, steps)
%
% Matrix multiplication version of trajectory sampling with a SPGP model.
%
% Inputs:
%   train_x: training inputs [N, d]
%   train_y: training targets [N, dim]
%   px: pseudo inputs [m, d]
%   hyp: log hyper parameters (b, c, sig) [d+2, dim]
%   start_states: starting states of trajectories to predict [n, dim]
%   polpar: !!!CURRENTLY LINEAR!!! policy parameters [dim, udim*n]
%   n: number of trajectories to predict [1]
%   steps: number of timesteps per trajectory to predict [1]
%
% Outputs:
%   trajs: the sampled trajectories
N = size(train_x, 1);
[m, d] = size(px);
[n, dim] = size(start_states);
udim = d - dim;

% Offilne computations
mean_factor = zeros(dim, m);
var_factor = zeros(m, m, dim);
sigf_sq_add_var = zeros(dim, n);
sign_sq_add_var = zeros(dim, n);
for i = 1:dim
    sigf_sq = exp(hyp(end-1, i));
    sign_sq = exp(hyp(end, i));
    Km = kern(px, px, hyp(:, i)) + 1e-6*eye(m);
    iKm = inv(Km);
    Kmn = kern(px, train_x, hyp(:, i));
    lamb = sigf_sq*ones(N, 1) - bsxfun(@dot, Kmn, Km\Kmn)';
    Qm = Km + Kmn*(bsxfun(@times, (1./(lamb + sign_sq))', Kmn))';
    mean_factor(i, :) = (Qm\(Kmn*(train_y(:, i)./(lamb + sign_sq))))';
    var_factor(:, :, i) = iKm - inv(Qm+1e-6*eye(m));
    
    sigf_sq_add_var(i, :) = sigf_sq*ones(1, n);
    sign_sq_add_var(i, :) = sign_sq*ones(1, n);
end

B = exp(hyp(1:end-2, :));
C = exp(hyp(end-1, :));
%px_ext = repmat(px, dim*n, 1); % [m, d] --> [m*n, d]
px_ext = repmat(px, n, 1); 
%C = repmat(C, m*n, 1); 
%C = reshape(C, 1, dim*m*n)'; 
su = repmat(start_states', udim, 1); su = reshape(su, dim, udim*n);
u = sum(bsxfun(@times, su, polpar), 1);

x = reshape([start_states, reshape(u, n, udim)]', 1, n*d);
s = start_states;

Broot = gpuArray(B.^.5);
B = gpuArray(B);
C = gpuArray(C);
px_ext = gpuArray(px_ext);
x = gpuArray(x);
xsave = parallel.gpu.GPUArray.zeros(steps, d*n);
xsave(1, :) = x;
mean_factor = gpuArray(mean_factor);
var_factor = gpuArray(var_factor);
sigf_sq_add_var = gpuArray(sigf_sq_add_var);
sign_sq_add_var = gpuArray(sign_sq_add_var);
s = gpuArray(s);
polpar = gpuArray(polpar);

for i = 2:steps
    % Online Computations
    x = reshape(x, d, n);
    x = repmat(x, m, 1); 
    x = reshape(x, d, m*n)'; 

	tmp1 = (px_ext - x);
	%tmp2 = tmp1.^2;

	for k = 1:dim
		tmp2 = arrayfun(@mysquare, tmp1, repmat(Broot(:, k)', m*n, 1), repmat(B(:, k)', m*n, 1));

		K = C(k)*exp(-.5*(tmp2*B(:, k)));
		K = reshape(K, m, n);

    	means(:, k) = (mean_factor(k, :)*K)';
    	vars(k, :) = sigf_sq_add_var(k, :) + sign_sq_add_var(k, :) - ...
        sum(bsxfun(@times, K, var_factor(:, :, k)*K), 1);
	end
    
    %vars = reshape(vars, dim, n)';
   
    ds = means + parallel.gpu.GPUArray.randn(n, dim).*(vars.^.5)';
    s = s+ds;
    %sflat = reshape(s', 1, dim*n);
	su = repmat(s', udim, 1); su = reshape(su, dim, udim*n);
    u = sum(bsxfun(@times, su, polpar), 1);
    x = reshape([s, reshape(u, udim, n)']', 1, n*d);
    xsave(i, :) = x;
    
end
xsave = gather(xsave);
