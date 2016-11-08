clear all, close all

r = {};
r{1} = load('814227_400ps_100p10_gold'); r{1} = r{1}.result;
r{2} = load('814228_400ps_40p4_gold'); r{2} = r{2}.result;
r{3} = load('814229_50ps_50p10_gold'); r{3} = r{3}.result;
r{4} = load('814230_100ps_50p10_gold'); r{4} = r{4}.result;
r{5} = load('814231_200ps_50p10_gold'); r{5} = r{5}.result;
r{6} = load('814232_300ps_50p10_gold'); r{6} = r{6}.result;
r{7} = load('814236_100ps_50p10_shortparallel'); r{7} = r{7}.result;
r{8} = load('814238_400ps_40p4_nogpucode_gold'); r{8} = r{8}.result;
r{9} = load('815341_200ps_100p10'); r{9} = r{9}.result;
r{10} = load('815342_200ps_100p10_shortpar_parfor'); r{10} = r{10}.result;
r{11} = load('815343_200ps_100p10_partforgpu'); r{11} = r{11}.result;
r{12} = load('815344_300ps_100p10_gpu'); r{12} = r{12}.result;
r{13} = load('test_wjvm_100ps_50p10'); r{13} = r{13}.result;
r{14} = load('test2_wjvm_100ps_50p10'); r{14} = r{14}.result;
r{15} = load('test3_wjvm_200ps_50p10'); r{15} = r{15}.result;
r{16} = load('test4_wjvm_300ps_50p10'); r{16} = r{16}.result;
r{17} = load('test5_wjvm_400ps_50p10'); r{17} = r{17}.result;
r{18} = load('test6_wjvm_200ps_50p10_siglik'); r{18} = r{18}.result;
r{19} = load('test7_wjvm_300ps_50p10_siglik'); r{19} = r{19}.result;
r{20} = load('test8_wjvm_500ps_50p10'); r{20} = r{20}.result;
r{21} = load('test9_wjvm_400ps_50p10_siglik'); r{21} = r{21}.result;
r{22} = load('test_wtraj_500ps50p10'); r{22} = r{22}.result;
r{23} = load('test_toeque5'); r{23} = r{23}.result;
r{24} = load('test_torque10'); r{24} = r{24}.result;
r{25} = load('test_torque15'); r{25} = r{25}.result;
r{26} = load('test_tq_q_noise1'); r{26} = r{26}.result;
r{27} = load('test_tq_q_noise2'); r{27} = r{27}.result;
r{28} = load('test_tq_q_noise3'); r{28} = r{28}.result;


% Time of one step of hyper parameter prediction vs number of pseudo inputs
timesh = [];
timest = [];
pseudos = [];
samples = [];
for i = 1:length(r)
	pseudos = [pseudos, r{i}.hyper_joint(1)*ones(1, length(r{i}.time.hyper))];
	samples = [samples, r{i}.new_sample_per_round:r{i}.new_sample_per_round:r{i}.total_samples];
	timesh = [timesh, r{i}.time.hyper./r{i}.hyper_joint(2)];
	timest = [timest, r{i}.time.pred./r{i}.total_samples];
end
figure, plot3(pseudos, samples, timesh, '*'), grid, xlabel('#pseudo inputs'), ylabel('#samples'), 
zlabel('calculation time [sec]'), title('Hyperparameter prediction one time step vs #pseudo inputs vs #samples')

% Prediction time of one trajectory vs number of pseudo inputs vs number of samples used
figure, plot3(pseudos, samples, timest, '*'), grid, xlabel('#pseudo inputs'), ylabel('#samples'), 
zlabel('prediction time [sec]'), title('prediction time of one trajectory vs #pseudo inputs vs #samples')

% Expected reward + variance + real reward for all samples
res = r{28}; 
for j = 1:uint8(res.total_samples/res.new_sample_per_round)
	figure(3)
	clf
	m1 = abs(res.rew(:, 1)-res.erew(:, 1, j))./abs(res.rew(:, 1));
	v1 = res.vrew(:, 1, j).^.5./abs(res.rew(:, 1));
	subplot(3,2,1)
	plot(log(m1), 'b'), hold on, ylabel('log torque penalty')
	plot(log(v1), 'r'),
	yLimits = get(gca,'YLim'); plot(double(j)*res.new_sample_per_round*ones(1, 2), yLimits, 'k--', 'LineWidth', 2)
	subplot(3,2,2)
	m2 = abs(res.rew(:, 2)-res.erew(:, 2, j))./abs(res.rew(:, 2));
	v2 = res.vrew(:, 2, j).^.5./abs(res.rew(:, 2));
	plot(log(m2), 'b'), hold on, ylabel('log torque constr')
	plot(log(v2), 'r'),
	yLimits = get(gca,'YLim'); plot(double(j)*res.new_sample_per_round*ones(1, 2), yLimits, 'k--', 'LineWidth', 2)
	subplot(3,2,3)
	m3 = abs(res.rew(:, 3)-res.erew(:, 3, j))./abs(res.rew(:, 3));
	v3 = res.vrew(:, 3, j).^.5./abs(res.rew(:, 3));
	plot(log(m3)), hold on, ylabel('log dist')
	plot(log(v3), 'r'),
	yLimits = get(gca,'YLim'); plot(double(j)*res.new_sample_per_round*ones(1, 2), yLimits, 'k--', 'LineWidth', 2)
	subplot(3,2,4)
	m4 = abs(res.rew(:, 4)-res.erew(:, 4, j))./abs(res.rew(:, 4));
	v4 = res.vrew(:, 4, j).^.5./abs(res.rew(:, 4));
	plot(log(m4), 'b'), hold on, ylabel('log q constr')
	plot(log(v4), 'r'),
	yLimits = get(gca,'YLim'); plot(double(j)*res.new_sample_per_round*ones(1, 2), yLimits, 'k--', 'LineWidth', 2)
	subplot(3,1,3)
	ms = abs(sum(res.rew(:, :), 2)-sum(res.erew(:, :, j), 2))./abs(sum(res.rew, 2));
	vs = sum(res.vrew(:, :, j), 2).^.5./abs(sum(res.rew, 2));
	plot(log(ms), 'b'), hold on, ylabel('log sum reward')
	plot(log(vs), 'r')
	yLimits = get(gca,'YLim'); plot(double(j)*res.new_sample_per_round*ones(1, 2), yLimits, 'k--', 'LineWidth', 2)
	title('expected squared distances between real and predicted individual and sum of rewards')

	figure(4)
	clf
	subplot(2, 1, 1)
	plot(log(res.joint_error_mean_var(:, 1, j)), 'b'), hold on,
	plot(log(res.joint_error_mean_var(:, 2, j).^.5), 'r')
	yLimits = get(gca,'YLim'); plot(double(j)*res.new_sample_per_round*ones(1, 2), yLimits, 'k--', 'LineWidth', 2)
	ylabel('log sq ball pos/vel error')
	subplot(2, 1, 2)
	plot(log(res.dist_error_mean_var(:, 1, j)), 'b'), hold on,
	plot(log(res.dist_error_mean_var(:, 2, j).^.5), 'r')
	yLimits = get(gca,'YLim'); plot(double(j)*res.new_sample_per_round*ones(1, 2), yLimits, 'k--', 'LineWidth', 2)
	ylabel('log sq distance error')

	figure(5)
    clf
    plotSampleTraj(6, res, j);
	
	pause
end
