function [xsave] = sampleTrajectoryCPULoopOverSamplesDimensions(train_x, train_y, px, hyp, ...
    start_states, polpar, steps)
%
% Matrix multiplication version of trajectory sampling with a SPGP model.
% !!! This version uses only CPU !!!
%
% Inputs:
%   train_x: training inputs [N, d]
%   train_y: training targets [N, dim]
%   px: pseudo inputs [m, d]
%   hyp: log hyper parameters (b, c, sig) [d+2, dim]
%   start_states: starting states of trajectories to predict [n, dim]
%   polpar: !!!CURRENTLY LINEAR!!! policy parameters [dim*n, udim*n]
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
bet_factor = zeros(dim, m*dim);
Lm_factor = zeros(dim*m, dim*m);
L_factor = zeros(dim*m, dim*m);
sigf_sq_add_var = zeros(1, dim*n);
sign_sq_add_var = zeros(1, dim*n);
for i = 1:dim
    [L, Lm, bet] = getLLmbet(train_x, train_y(:, i), px, hyp(:, i));
    ix = ((i-1)*m+1):(i*m);
    Lm_factor(ix, ix) = Lm;
    L_factor(ix, ix) = L;
    bet_factor(i, ix) = bet';
    sigf_sq = exp(hyp(end-1, i));
    sign_sq = exp(hyp(end, i));
    ix = i:dim:(dim*n);
    sigf_sq_add_var(ix) = sigf_sq*ones(1, n);
    sign_sq_add_var(ix) = sign_sq*ones(1, n);
end
B = exp(hyp(1:end-2, :));
C = exp(hyp(end-1, :));
px_ext = repmat(px, dim*n, 1); % [m, d] --> [m*n*dim, d]
B = repmat(B, m, 1); %[d, dim] --> [m*d, dim]
B = reshape(B, d, m*dim)'; 
B = repmat(B, n, 1); 
C = repmat(C, m, 1); % [1, dim] --> [m, dim]
C = reshape(C, m*dim, 1); % [m, dim] -> [m*dim, 1]
C = repmat(C, n, 1); % [n*m*dim, 1]
sflat = reshape(start_states', 1, dim*n);
u = sflat*polpar;
x = reshape([start_states, reshape(u, n, udim)]', 1, n*d);
s = start_states;
xsave = zeros(steps, d*n);
xsave(1, :) = x;

for i = 2:steps
    % Online Computations
    x = repmat(x, dim*m, 1); % [1, n*d] --> [dim*m, n*d]
    x = reshape(x', d, m*n*dim)'; % [dim*m, n*d] --> [m*n*dim, d]

    clear K
    K = C.*exp(-.5*sum(B.*(px_ext - x).^2, 2));
    K = reshape(K, dim*m, n);
    for j = 1:n
        for k = 1:dim
            ixd = (k-1)*m+1:k*m;
            lst(ixd, j) = L_factor(ixd, ixd)\K(ixd, j);
            lmst(ixd, j) = Lm_factor(ixd, ixd)\lst(ixd, j);
            means(j, k) = (bet_factor(k, ixd)*lmst(ixd, j))';
            ix = (j-1)*dim+1:(j-1)*dim+k;
            vars(ix) = - sum(lst(ixd, j).^2, 1) ... 
           + exp(hyp(end, k)).*sum(lmst(:, j).^2, 1) ...
           + exp(hyp(end-1, k)) + exp(hyp(end, k));
        end
    end
    vars = reshape(vars, dim, n)';
   
    rng(i)
    ds = means + randn(n, dim).*(vars.^.5);
    s = s+ds;
    sflat = reshape(s', 1, dim*n);
    u = sflat*polpar;
    x = reshape([s, reshape(u, udim, n)']', 1, n*d);
    xsave(i, :) = x;
    
end


