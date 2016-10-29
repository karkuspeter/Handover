clear all, close all

addpath('../gpml')
addpath('../SPGP_dist')
startup % of GPML

sig_noise = [.1, .2];
K = -[3.1, 1.7];
x0 = [0, 1];
dt = .01;
diffx_int2 = @(t, x, k) [x(2), x(:)'*k(:)]';
diffx_loc = @(t, x) diffx_int2(t, x, K);

% Generate observations
tm = [0:dt:10];
xsave = [];
xsaveex = [];
usave = [];
notrajs = 10;
for e = 1:notrajs
	x = 2*(rand(1, 2) -.5)';
	y = zeros(length(tm), 2);
	for t = 2:length(tm)
		x = x + diffx_loc(0, x)*dt + randn(2, 1).*sig_noise'*dt;
		y(t, :) = x';
	end
	xsave = [xsave; y];
	xsaveex = [xsaveex y];
	usave = [usave; y*K'];
end

% Prepare training data
trainx = [xsave, usave];
train_y = [];
train_x = [];
for i = 1:notrajs
	train_y = [train_y; diff(trainx((i-1)*length(tm)+2:i*length(tm), 1:2))];
	train_x = [train_x; trainx( (i-1)*length(tm)+2:i*length(tm)-1, :)];
end

% Optimize hyper/pseudo
m = 100; 
[hyp, px] = hyper_pseudo_commonInputs(train_x, train_y, m, 150); 

% transform hyp log(w, sigf, sign) to log(b, sigf^2, sign^2)
for i = 1:2
    w = exp(hyp(1:end-2, i));
    sigf = exp(hyp(end-1, i));
    sign = exp(hyp(end, i));
    hyp(:, i) = log([w.^-2; sigf^2; sign^2]);
end

%Precalculate parameters
m = size(px, 1);
L = zeros(m, m, 2); Lm = zeros(m, m, 2); bet = zeros(m, 2);
for i = 1:2 % state dimension
    [d1, d2, d3] = getLLmbet(train_x, train_y(:, i), px, hyp(:, i));
    L(:, :, i) = d1;
    Lm(:, :, i) = d2;
    bet(:, i) = d3;
end

x0 = [1 0];
K = -[3.1, 1.7];
notrajs = 5; 
len = length(tm)-2;
trajs = zeros(len, 3, notrajs);
figure,
try
	matlabpool
catch
	disp('matlabpool already open')
end
tic
parfor i = 1:notrajs
   traj = sampleTrajectory(x0, L, Lm, bet, px, hyp, K', len);
   trajs(:, :, i) = traj;
   plot(traj), hold on
end
time_cpu = toc;

disp(['CPU parallel: ',num2str(time_cpu)])

traj_result.cpu_parallel = trajs;

notrajs = 10;
K = -[3.1, 1.7];
polpar = K;
K = zeros(2*notrajs, notrajs);
K(1:2:end, 1:end) = diag(polpar(1)*ones(notrajs, 1));
K(2:2:end, 1:end) = diag(polpar(2)*ones(notrajs, 1));
start_states = repmat([0, 1], notrajs, 1);

tic
trajs = sampleTrajectoryCPU(train_x, train_y, px, hyp, start_states, K, len);
time_cpu_matrix = toc;
disp(['CPU matrix: ',num2str(time_cpu_matrix)])

% traj_result.cpu_matrix = trajs;

tic
trajs = sampleTrajectoryCPULoopOverDimensions(train_x, train_y, px, hyp, start_states, K, len);
time_cpu_matrix_loop = toc;
disp(['CPU matrix dimloop: ',num2str(time_cpu_matrix_loop)])
% traj_result.cpu_matrix = trajs;

tic
trajs = sampleTrajectoryCPU2(train_x, train_y, px, hyp, start_states, K, len);
time_cpu_matrix_ver2 = toc;
disp(['CPU matrix ver2: ',num2str(time_cpu_matrix_ver2)])
% traj_result.cpu_matrix = trajs;

% tic
% trajs = sampleTrajectoryGPU(train_x, train_y, px, hyp, start_states, K, len);
% time_gpu_matrix = toc;
% disp(['GPU matrix: ',num2str(time_gpu_matrix)])
% traj_result.gpu_matrix = trajs;
