clear all, close all

addpath('../gpml')
addpath('../SPGP_dist')
startup % of GPML

trajs = load('traj5to15.mat');
trajs = trajs.traj;
len = trajs.len;

xtrain = [];
ytrain = [];

q = trajs.q;
u = trajs.u;
ref = trajs.ref;

notrajs = size(q, 1)/len;

for i = 1:notrajs

	ix = ((i-1)*len+1) : (i*len -1);
	transq = q(ix, :);
	% add 2pi to the last joint angle to come over the gap during prediction
	ixtrans = find(transq(:, 7) < 0);
	transq(ixtrans, 7) = transq(ixtrans, 7) + 2*pi;  
	xtrain = [xtrain; transq, u(ix, :)];
	ixx = ((i-1)*len+1) : (i*len);
	ytrain = [ytrain; diff(q(ixx, :))];
end
ytrain(:, 1:2:end) = mpi2pi(ytrain(:, 1:2:end));

% Truncate database
% G = 2000;
% ixok = randperm(size(ytrain, 1));
% xtrain = xtrain(1:len, :);
% ytrain = ytrain(1:len, :);

if 1
	[hyp, px] = hyper_pseudo_commonInputs(xtrain, ytrain, 50, 300); 
else
model = load('model.mat');
	hyp = model.model.hyp;
	px = model.model.px;
end

% Extra parameters from simulation
PGains = trajs.props.environment.PGains;
DGains = trajs.props.environment.DGains;

% Prediction parameters
K = zeros(4, 8);
pg = -diag(PGains);
dg = -diag(DGains);
K(1, 1:2) = [pg(1), dg(1)];
K(2, 3:4) = [pg(2), dg(2)];
K(3, 5:6) = [pg(3), dg(3)];
K(4, 7:8) = [pg(4), dg(4)];

control_input = @(x, k, r) k*(x(:)-r(:));
iscov = @(S, k) k*S;
iicov = @(S, k) k*S*k';


% Extracting Reference Trajectory
dt = trajs.props.dt;
% !!! transform angle data to avoid big jumps in velocity !!!
reftr = ref(1:len, :); reftr(reftr(:, 4) < 0, 4) = reftr(reftr(:, 4) < 0, 4) + 2*pi;
[pref, dref] = extract_velocity(ref(1:len, :), dt);
ref = zeros(size(pref, 1), 8);
ref(:, 1:2:end) = pref;
ref(:, 2:2:end) = dref;
m = ref(1, :)';
s = 1e-6*eye(8);

[ms, covs, sigs] = trajGpPredPseudo(hyp, xtrain, px, ytrain, m, s, ref, control_input, iscov, iicov, K);



