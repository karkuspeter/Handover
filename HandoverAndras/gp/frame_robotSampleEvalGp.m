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
N = 6; % Generate this much trajectories
M = 2; % Inlcude this much trajectories per round
dt_projectile = .1;
pseudos = 10;
pseudos_release = M;
pseudos_projectile = M;
hyp_opt_steps = 50;
hyp_opt_steps_release = 20;
hyp_opt_steps_projectile = 20;
properties.environment.noise = 10*ones(1, 4);
% ------------------

%% Generating trajectories
pseudos_release = min(pseudos_release, N);
pseudos_projectile = min(pseudos_projectile, N);
times = {};
tic
x = [];
y = [];
ref = [];
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
keyboard
len = size(q, 2) - 1;
y(:, 1:2:end) = mpi2pi(y(:, 1:2:end));
times.traj_gen = toc;
times.hyper = [];
times.pred = [];

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
for i = M:M:N

	hyper_start = toc;
    if 1
		% Always start hyperparameter optimization from beginning
        [hyp, px, hyper_traj] = hyper_pseudo_commonInputs(x(1:(i*len), :), y(1:(i*len), :), pseudos, hyp_opt_steps); 
		[hyp_release, px_release] = hyper_pseudo_commonInputs(xrelease(1:i, :), yrelease(1:i, :), i, hyp_opt_steps_release);
		[hyp_projectile, px_projectile] = hyper_pseudo_commonInputs(xprojectile(1:i, :), yprojectile(1:i, :), i, hyp_opt_steps_projectile);
    else
		% Reuse previously optimized hyperparameters 
        [hyp, px] = hyper_pseudo_commonInputs(x(1:(i*len), :), y(1:(i*len), :), 300, 50, hyp, px); 
    end
    hyper_traj_save(:, :, uint8(i/M)) = hyper_traj;
    
	hyper_end = toc;
	times.hyper = [times.hyper, hyper_end-hyper_start];

	gpmodel_release.inputs = xrelease(1:i, :);
	gpmodel_release.target = yrelease(1:i, :);
	gpmodel_release.hyp = hyp_release;
	gpmodel_release.induce = px_release;
	
	gpmodel_proj.inputs = xprojectile(1:i, :);
	gpmodel_proj.target = yprojectile(1:i, :);
	gpmodel_proj.hyp = hyp_projectile;
	gpmodel_proj.induce = px_projectile;

	pred_start = toc;
	% Test on all samples
    dummy1 = zeros(N, 2);
    dummy2 = zeros(N, 2);
    erew_loc = zeros(N, 4); vrew_loc = zeros(N, 4);
    sq_track_error_mean_loc = zeros(N, 12);
    sq_track_error_var_loc = zeros(N, 12);
    
    traj_mean_save_loc = zeros(len, 12, N);
    traj_var_save_loc = zeros(len, 12, N);

	parfor j = 1:N

		ix = ((j-1)*len+1):(j*len);
		ref_new = ref(ix, :);
		m = ref_new(1, :)';
		s = 1e-6*eye(8);

		% Traj prediction
		% ... joint trajectory
		[ms, covs, sigs] = trajGpPredPseudo(hyp, x(1:(i*len), :), px, ...
			y(1:(i*len), :), m, s, ref_new, control_input, iscov, iicov, K);
        traj_mean_save_loc(:, :, j) = ms;
        traj_var_save_loc(:, :, j) = sigs;
		% ... init projectile pos/vel
		[mproj, covproj, sigproj] = gp1(gpmodel_release, ms(ix_release(j), 1:8)', covs(1:8, 1:8, ix_release(j)));
		covproj = covproj + diag(exp(2*hyp_release(end, :)));
		% ... projectile distance
		[mdist, covdist, sigdist] = gp1(gpmodel_proj, mproj, covproj);
		covdist = covdist + diag(exp(2*hyp_projectile(end, :)));


		% Expected reward
		% ... from torque, torque constraint, joint constraint
		[mu_rew_torque, var_rew_torque] = getQuadraticCost(ms(:, 9:12), covs(9:12, 9:12, :), -5e-6*dt*eye(4));
		[mu_rew_torqueViolation, var_rew_torqueViolation] = getConstraintViolationCost(ms(:, 9:12), covs(9:12, 9:12, :), model.minU', model.maxU', -1e-1*dt*eye(4));
		[mu_rew_jointViolation, var_rew_jointViolation] = getConstraintViolationCost(ms(:, 1:8), covs(1:8, 1:8, :), model.minQ', model.maxQ', -1e5*dt*eye(8));
		% ... from hitting target
		[mu_rew_distance, var_rew_distance] = getQuadraticCost(mdist - states(j, 1), covdist, -100);
		
		erew_loc(j, :) = [mu_rew_torque mu_rew_torqueViolation mu_rew_distance mu_rew_jointViolation];
		vrew_loc(j, :) = [var_rew_torque var_rew_torqueViolation var_rew_distance var_rew_jointViolation];

		% Square Joint Tracking Error
		for k = 1:12
			[mut, vart] = getQuadraticCost(ms(:, k)-x(ix, k), covs(k, k, :));
			sq_track_error_mean_loc(j, k) = mut;
			sq_track_error_var_loc(j, k) = vart;
		end
		% Square joint -> starting pos/vel error
		[mut, vart] = getQuadraticCost(mproj' - yrelease(j, :), covproj);
        dummy1(j, :) = [mut, vart];
% 		sq_joint_error_mean_var_loc(j, 1) = mut;
% 		sq_joint_error_mean_var_loc(j, 2) = vart;
		% Starting pos/vel error --> distance error
		[mut, vart] = getQuadraticCost(mdist - yprojectile(j), covdist);
		dummy2(j, :) = [mut vart];
% 		sq_dist_error_mean_var_locvar(j, 2) = vart;
% 		sq_dist_error_mean_var_locvar(j, 1) = mut;
    
    end
    % Merging with big database 
    tmpi = uint8(i/M);
    erew(:, :, tmpi) = erew_loc;
    vrew(:, :, tmpi) = vrew_loc;
    sq_track_error_mean(:, :, tmpi) = sq_track_error_mean_loc;
    sq_track_error_var(:, :, tmpi) = sq_track_error_var_loc;
    sq_joint_error_mean_var(:, :, tmpi) = dummy1;
    sq_dist_error_mean_var(:, :, tmpi) = dummy2;
    
    for j = 1:N
        ix = ((j-1)*len+1):(j*len); 
        traj_mean_save(ix, :, tmpi) = traj_mean_save_loc(:, :, j);
        traj_var_save(ix, :, tmpi) = traj_var_save_loc(:, :, j);
    end
  
    
	pred_end = toc;
	times.pred = [times.pred, pred_end-pred_start];
    
end
result.new_sample_per_round = M;
result.total_samples = N;
result.rew = rew;
result.erew = erew;
result.vrew = vrew;
result.track_error_mean = sq_track_error_mean;
result.track_error_var = sq_track_error_var;
result.time = times;
result.joint_error_mean_var = sq_joint_error_mean_var;
result.dist_error_mean_var = sq_dist_error_mean_var;
result.hyper_joint = [pseudos, hyp_opt_steps];
result.hyper_release = [pseudos_release, hyp_opt_steps_release];
result.hyper_projectile = [pseudos_projectile, hyp_opt_steps_projectile];
result.x = x;
result.len = len;
result.hyper_traj = hyper_traj_save;
result.traj_mean = traj_mean_save;
result.traj_std = traj_var_save;
save(name_str, 'result');
try 
	matlabpool close
catch
	disp('no open matlabpools')
end
