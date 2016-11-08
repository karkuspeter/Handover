function [xsave] = sampleTrajectoryCPU2(px, hyp, mean_factor, ...
    var_factor, start_states, polpar, s, steps)
%
% Matrix multiplication version of trajectory sampling with a SPGP model.
% !!! This version uses only CPU !!!
%
% Inputs:
%   px: pseudo inputs [m, d]
%   hyp: log hyper parameters (b, c, sig) [d+2, dim]
%	mean_factor: parameters for mean prediction [dim, m] (see sampleTrajectoryPrecalculations.m)
%	var_factor: parameters for mean prediction [m, m, dim] (see sampleTrajectoryPrecalculations.m)
%   start_states: starting states of trajectories to predict [samples, dim]
%   polpar: the linear policy parameters [udim, dim*samples], s.t. u = polpar*(ref - q)
%   s: number of trajectories to predict per sample [1]
%   steps: number of timesteps per trajectory to predict [1]
%		    --- OR ----
%		  a reference trajectory for states [steps, samples*dim] 
% Outputs:
%   trajs: the sampled trajectories

[m, d] = size(px);
[samples, dim] = size(start_states);
udim = d - dim;

n = s*samples;
start_states = reshape(repmat(start_states', s, 1), dim, s*samples)';
polpar_new = zeros(n*udim, dim);
for i = 1:samples
	polpar_new((i-1)*s*udim+1:i*s*udim, :) = repmat(polpar(:, (i-1)*dim+1:i*dim), s, 1);
end
polpar = polpar_new';

if size(steps, 1) > 1
	tracking = 1;
	ref = steps;
	steps = size(ref, 1);
else
	tracking = 0;
end

% Offilne computations
Mfs = zeros(m, 1, dim);
for i = 1:dim
    sigf_sq = exp(hyp(end-1, i));
    sign_sq = exp(hyp(end, i));
    sigf_sq_add_var(i, :) = sigf_sq*ones(1, n);
    sign_sq_add_var(i, :) = sign_sq*ones(1, n);
	Mfs(:, :, i) = mean_factor(i, :)';
end
B = exp(hyp(1:end-2, :));
C = exp(hyp(end-1, :));
px_ext = repmat(px, n, 1); 

su = reshape(repmat(start_states', udim, 1), dim, udim*n);
if ~tracking
	u = sum(bsxfun(@times, su, polpar), 1);
else
	ref_loc = reshape(repmat(reshape(ref(1, :), dim, samples), s, 1), dim, n);
	ref_loc = reshape(repmat(ref_loc, udim, 1), dim, udim*n);
	u = sum(bsxfun(@times, su - ref_loc, polpar), 1);
end

x = reshape([start_states, reshape(u, n, udim)]', 1, n*d);
ss = start_states;

xsave = zeros(steps, d*n);
xsave(1, :) = x;

for i = 2:steps
    % Online Computations
    x = reshape(x, d, n);
    x = repmat(x, m, 1); 
    x = reshape(x, d, m*n)'; 

	tmp1 = (px_ext - x);
	tmp2 = tmp1.^2;

	if 0 % Paralellized version
		Kd = bsxfun(@times, exp(-.5*tmp2*B), C);
		Kd = reshape(Kd, [m, n, dim]);
		means = squeeze(sum(bsxfun(@times, Kd, Mfs), 1));
		for k = 1:dim
			vars(k, :) = sigf_sq_add_var(k, :) + sign_sq_add_var(k, :) - ...
			sum(bsxfun(@times, Kd(:, :, k), var_factor(:, :, k)*Kd(:, :, k)), 1);
		end
	else
		for k = 1:dim

			K = C(k)*exp(-.5*(tmp2*B(:, k)));
			K = reshape(K, m, n);

			means(:, k) = (mean_factor(k, :)*K)';
			vars(k, :) = sigf_sq_add_var(k, :) + sign_sq_add_var(k, :) - ...
			sum(bsxfun(@times, K, var_factor(:, :, k)*K), 1);
		end
	end
		
   
    ds = means + randn(n, dim).*(vars.^.5)';
    ss = ss+ds;
	su = repmat(ss', udim, 1); su = reshape(su, dim, udim*n);
	if ~tracking
		u = sum(bsxfun(@times, su, polpar), 1);
	else
		ref_loc = reshape(repmat(reshape(ref(i, :), dim, samples), s, 1), dim, n);
		ref_loc = reshape(repmat(ref_loc, udim, 1), dim, udim*n);
		u = sum(bsxfun(@times, ref_loc - su, polpar), 1);
	end
    x = reshape([ss, reshape(u, udim, n)']', 1, n*d);
    xsave(i, :) = x;
 
end
