%% Init
initQuadLink; 
dt = model.dt;
% Savename generation
curr_clock = clock;
dummy = curr_clock(6)*10000;
dummy = num2str(dummy);
name_str = [date,'-',num2str(curr_clock(4)),'-',num2str(curr_clock(5)),'-',dummy(1:6),'.mat'];
try
	matlabpool
catch
	disp('matlabpool already open')
end
%% --- PARAMETERS -----
N = 200; % Generate this much trajectories
M = 200;
dt_projectile = .1;
pseudos = 400;
pseudos_release = N;
pseudos_projectile = N;
hyp_opt_steps = 150;
hyp_opt_steps_release = 200;
hyp_opt_steps_projectile = 200;
tqnoise = 3;
state_noise = [.02, .2, .02, .2, .02, .2, .02, .2]./3;
% ------------------
%% Generating trajectories
properties.environment.noise=tqnoise*ones(1, 4);
properties.environment.state_noise = state_noise;
pseudos_release = min(pseudos_release, N);
pseudos_projectile = min(pseudos_projectile, N);
times = {};
tic
x = [];
y = [];
ref = [];
ref_learn = []; % maybe it's not necessary but I don't want to change ref now.
rew = [];
xrelease = [];
yrelease = [];
xprojectile = [];
yprojectile = [];
ix_release = [];
states = zeros(N, 2);
for i = 1:N
    rng(i) % set seed
    state = rand(1, 2).*[10, 0] + [5, 0];
	states(i, :) = state;
    sample = model.sampleFunc(model, state, 1);
    [reward, q, u, yGen, yGend, traj_extra] = properties.rewardFunction([], properties, sample);
	rewDist = traj_extra.rewardDist;
    transq = q';
	% add 2pi to the last joint angle to come over the gap during prediction
	ixtrans = find(transq(1:end-1, 7) < 0);
	transq(ixtrans, 7) = transq(ixtrans, 7) + 2*pi;  
	x = [x; transq(1:end-1, :), u(:, 1:end-1)'];
	y = [y; diff(q')];
    ref_loc = zeros(size(yGen, 1), 8);
    ref_loc(:, 1:2:end) = yGen;
    ref_loc(:, 2:2:end) = yGend;
    ref = [ref; ref_loc(1:end-1, :)];
	ref_learn = [ref_learn ref_loc(1:end-1, :)];
	rew = [rew; rewDist];
	qrelease = transq(floor(sample(end)/dt), :);
	xrelease = [xrelease; qrelease];
	[J, startpos] = getJacobian(qrelease(1:2:end), properties.environment, 4);
	yrelease = [yrelease; startpos', qrelease(2:2:end)*J'];
	ix_release = [ix_release, traj_extra.releaseIndex];
	%[xprojectile_curr, yprojectile_curr] = getProjectileTrajectory(yrelease(end, :), dt_projectile);
	%xprojectile = [xprojectile; xprojectile_curr];
	%yprojectile = [yprojectile; yprojectile_curr];
	%len_projectile = [len_projectile, size(xprojectile, 1)];
	d = getProjectileLandingDistance(yrelease(end, :));
	xprojectile = [xprojectile; yrelease(end, :)];
	yprojectile = [yprojectile; d];
end
len = size(q, 2) - 1;
y(:, 1:2:end) = mpi2pi(y(:, 1:2:end));
times.traj_gen = toc;
%res.x = x;
%res.y = y;
%save res.mat res;
%keyboard

tmp = uint8(N/M);
erew = zeros(N, 4, tmp); vrew = zeros(N, 4, tmp);
sq_track_error_mean = zeros(N, 12, tmp);
sq_track_error_var = zeros(N, 12, tmp);
sq_joint_error_mean_var = zeros(N, 2, tmp); 
sq_dist_error_mean_var = zeros(N, 2, tmp); 

traj_mean_save = zeros(len*N, 12, tmp);
traj_var_save = zeros(len*N, 12, tmp);
hyper_traj_save = zeros(hyp_opt_steps, pseudos*12 + 8*14 + 1, tmp);


%% Hyperparameter optimization and prediction
hyper_start = toc;
if 0
		% Always start hyperparameter optimization from beginning
        [hyp, px] = hyper_pseudo_commonInputs(x, y, pseudos, hyp_opt_steps); 
		[hyp_release, px_release] = hyper_pseudo_commonInputs(xrelease, yrelease, N, hyp_opt_steps_release);
		[hyp_projectile, px_projectile] = hyper_pseudo_commonInputs(xprojectile, yprojectile, N, hyp_opt_steps_projectile);
		robot_hyperpx.hyp = hyp;
		robot_hyperpx.px = px;
		robot_hyperpx.hyp_release = hyp_release;
		robot_hyperpx.px_release = px_release;
		robot_hyperpx.hyp_projectile = hyp_projectile;
		robot_hyperpx.px_projectile = px_projectile;
		save robot_hyperpx.mat robot_hyperpx
else
	dummy = load('robot_hyperpx');
	robot_hyperpx = dummy.robot_hyperpx;
	hyp = robot_hyperpx.hyp;
	px = robot_hyperpx.px;
	hyp_release = robot_hyperpx.hyp_release;
	px_release = robot_hyperpx.px_release;
	hyp_projectile = robot_hyperpx.hyp_projectile;
	px_projectile = robot_hyperpx.px_projectile;
end


% transform hyp log(w, sigf, sign) to log(b, sigf^2, sign^2)
for i = 1:size(hyp, 2)
    w = exp(hyp(1:end-2, i));
    sigf = exp(hyp(end-1, i));
    sign = exp(hyp(end, i));
    hyp(:, i) = log([w.^-2; sigf^2; sign^2]);
end

% -------------- SIMULATION PARAMETERS -----------------
notrajs = 200; % reproduce this many trajectories...
s = 50;      % ...this many times
start_states = repmat([0, 0, 0, 0, 0, 0, pi, 0], notrajs, 1);
polpar = zeros(4, 8);
polpar(:, 1:2:end) = properties.environment.PGains;
polpar(:, 2:2:end) = properties.environment.DGains;
polpar  = repmat(polpar, 1, notrajs);
ix = 1:notrajs*8;
ref = ref_learn(:, ix);

tic
[mean_factor, var_factor] = sampleTrajectoryPrecalculations(x, y, px, hyp);
times.precalculations = toc;
reset(gpuDevice)

% ------------------ CPU2 ---------------------------
if 1
%profile on
	workers = matlabpool('size');
	tic
	trajs = sampleTrajectoryCPU2(px, hyp, mean_factor, var_factor, start_states, polpar, s, ref);
	time_cpu_matrix = toc;
	times.cpu = time_cpu_matrix;
	disp(['CPU2: ',num2str(time_cpu_matrix)])
%profsave(profile('info'), 'profsave_cpu2');
%	traj_result.cpu2 = trajs;
end


% ------------------ CPU2 parallel ---------------------------
if 1
%profile on
	workers = matlabpool('size');
	tic
	spmd(workers)
		trajs = sampleTrajectoryCPU2(px, hyp, mean_factor, var_factor, start_states, polpar, ceil(s/workers), ref);
	end
	time_cpu_matrix = toc;
	times.cpu_parallel = time_cpu_matrix;
	disp(['CPU2 parallel: ',num2str(time_cpu_matrix)])
%profsave(profile('info'), 'profsave_cpu2');
%	traj_result.cpu2 = trajs;
end

% ------------------- GPU3 ------------------
if 1
%profile on
	tic
	trajs = sampleTrajectoryGPU3(px, hyp, mean_factor, var_factor, start_states, polpar, s, ref);
	time_gpu_matrix = toc;
	times.gpu_efficient = time_gpu_matrix;
	disp(['GPU3: ',num2str(time_gpu_matrix)])
%profsave(profile('info'), 'profsave_gpu3');
%	traj_result.cpu2 = trajs;
end

% ------------------- GPU3 single ------------------
if 0
%profile on
	tic
	trajs = sampleTrajectoryGPU3single(px, hyp, mean_factor, var_factor, start_states, polpar, s, ref);
	time_gpu_matrix = toc;
	times.gpu_efficient = time_gpu_matrix;
	disp(['GPU3 single: ',num2str(time_gpu_matrix)])
%profsave(profile('info'), 'profsave_gpu3');
%	traj_result.cpu2 = trajs;
end

% ------------------- GPU4 ------------------
if 0
%profile on
	tic
	trajs = sampleTrajectoryGPU4(px, hyp, mean_factor, var_factor, start_states, polpar, s, ref);
	time_gpu_matrix = toc;
	times.gpu_own_square = time_gpu_matrix;
	disp(['GPU4: ',num2str(time_gpu_matrix)])
%profsave(profile('info'), 'profsave_gpu4');
%	traj_result.cpu2 = trajs;
end


% ------------------- GPU3 paralell in dim ------------------
if 0
profile on
	tic
	trajs = sampleTrajectoryGPU3parallelDim(px, hyp,mean_factor, var_factor, start_states, polpar, s, ref);
	time_gpu_matrix = toc;
	times.gpu_par_dim = time_gpu_matrix;
	disp(['GPU3 par dim: ',num2str(time_gpu_matrix)])
profsave(profile('info'), 'profsave_gpu3paralleldim');
%	traj_result.cpu2 = trajs;
end

% ------------------ PARALLELL GPU3 -------------------------------
if 1

noGPUdevice = gpuDeviceCount;
try
	matlabpool
catch
	disp('matlabpool already open')
end
tic
spmd(noGPUdevice)
	i = labindex;
	gpuDevice(i);
	trajs = sampleTrajectoryGPU3(px, hyp, mean_factor, var_factor, start_states, polpar, round(s/noGPUdevice), ref);
end
time_gpu_parallel = toc;
times.gpu_efficient_multigpu = time_gpu_parallel;
disp(['GPU3 parallel: ',num2str(time_gpu_parallel)])
% Getting data from composite
%trajs_save = [];
%for i = 1:noGPUdevice
%	trajs_save = [trajs_save; trajs{i}];
%end

% traj_result.gpu2 = trajs;
% save traj_result.mat traj_result
end
times
