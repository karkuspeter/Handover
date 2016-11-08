function [trajs] = sampleTrajectoryCPU(train_x, train_y, px, hyp, n, steps);
%
% Matrix multiplication version of trajectory sampling with a SPGP model.
% !!! This version uses only CPU !!!
%
% Inputs:
%   train_x: training inputs [N, d]
%   train_y: training targets [N, dim]
%   px: pseudo inputs [m, d]
%   hyp: log hyper parameters (b, c, sig) [d+2, dim]
%   n: number of trajectories to predict [1]
%   steps: number of timesteps per trajectory to predict [1]
%
% Outputs:
%   trajs: the sampled trajectories

% Offilne computations
mean_factor = zeros(m, dim);
Qm_factor = zeros(m, dim*m);
Km_factor = zeros(m, dim*m);
for i = 1:dim
    sigf_sq = exp(hyp(end-1, i));
    sign_sq = exp(hyp(end, i));
    Km = kern(px, px, hyp(:, i)); % [m, m]
    Knm = kern(px, train_x, hyp(:, i)); % [m, n]
    lamb = sigf_sq*ones(n, 1) - diag(Knm' * Km\Knm); % Km\Knm = inv(Km)*Knm
    Gamma = diag(lamb);
    Qm = Km + Knm/(Gamma + sig_sq*eye(n))*Knm';
    dummy = Knm*(Gamma + sig_sq*eye(n))*train_y(:, i);
    mean_factor(:, i) = Qm\dummy;
    Km_factor(:, ((i-1)*m+1):(i*m)) = Km;
    Qm_cactor(:, ((i-1)*m+1):(i*m)) = Qm;
end
B = exp(hyp(1:end-2, :));
C = exp(hyp(end-1, :));
px_ext = repmat(px, dim*n, 1); % [m, d] --> [m*n*dim, d]
B = repmat(B, m, 1); %[d, dim] --> [m*d, dim]
B = reshape(B, m*d*dim, 1); % [m*d, dim] --> [m*d*dim, 1]
B = repmat(B, 1, d); % [m*d*dim, 1] --> [m*d*dim, d]
C = repmat(C, m, 1); % [1, dim] --> [m, dim]
C = reshape(C, m*dim, 1); % [m, dim] -> [m*dim, 1]
C = repmat(C, n, 1); % [n*m*dim, 1]

for i = 1:steps
    % Online Computations
    x = repmat(x, dim*m, 1); % [1, n*d] --> [dim*m, n*d]
    x = reshape(x, m*n*dim, d); % [dim*m, n*d] --> [m*n*dim, d]

    K = C.*exp(sum(-.5*B.*(px_ext - x_ext).^2, 2));
end

