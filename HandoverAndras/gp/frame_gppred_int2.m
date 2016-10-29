clear all, close all

K = -[3, 1];
x0 = [0, 1];
dt = .01;
dtgp = .01;
tm = [0:dt:10];

gpstep = dtgp/dt;

diffx_int2 = @(t, x, k) [x(2), x(:)'*k(:)]';

% Prediction parameters
control_input = @(x, k) x(:)'*k(:);
iscov = @(S, K) K(:)'*S;
iicov = @(S, K) K(:)'*S*K(:);
limit_state = @(x) x;

allData = [];
figure,
for e = 1:10

	x0 = (rand(1, 2)-.5)*2;

	tm = [0:dt:10];
	K = -[2; 3] + 2*(rand(2, 1)-.5).*[1; 2];
	diffx_loc = @(t, x) diffx_int2(t, x, K);
	%[t,y] = ode45(diffx_loc, tm, x0);
	x = x0';
	y = zeros(length(tm), 2);
	y(1, :) = x';
	for t = 2:length(tm)
		x = x + diffx_loc(0, x)*dt + randn(2, 1).*[.1; .1]*dt;%*(t*dt)^.5;
		y(t, :) = x';
	end

	s = y;
	a = y*K(:);

	train_x = [s(1:gpstep:end-1, :), a(1:gpstep:end-1)];
	train_y = [diff(s(1:gpstep:end, :))];

	allData = [allData; train_x];
	
	% Sparsification
	for i = 1:size(train_x, 1)
		dnew = [train_x(i, :), train_y(i, :)];
		[S, update, ix_update] = incr_online_sparsification(dnew, thresh, N, kernel_pars, S, 0);
	end

	% Hyperparameter learning
	options = optimset('GradObj', 'on');
	options = optimset(options, 'Display', 'off');
	for i = 1:2
		[hyp_opt, fval] = fminunc(@(hyp) hyper_optim(hyp, S.D(:, 1:3), S.D(:, 3+i), 0), Hyp(i, :), options);
		Hyp(i, :) = hyp_opt;
	end

	% Model validation
	[logl, ll] = validateGpModel(Hyp, S.D(:, 1:3), S.D(:, 4:5), 0);

	% Prediction with same starting sates and gain
	[ms, covs, sigs] = trajgpPred(Hyp, S.D(:, 1:3), S.D(:, 4:5), x0', diag([1e-6, 1e-6]), (10/dt/gpstep)+1, control_input, iscov, iicov, K(:)', limit_state);

	% Prediction with same starting states and other gain
%	tm = [0:dt:10];
%	K1 = -[2; 3] + 2*(rand(2, 1)-.5).*[1; 2];
%	diffx_loc = @(t, x) diffx_int2(t, x, K1);
%	%[t,y] = ode45(diffx_loc, tm, x0);
%	x1 = x0';
%	y1 = zeros(length(tm), 2);
%	y1(1, :) = x1';
%	for t = 2:length(tm)
%		x1 = x1 + diffx_loc(0, x1)*dt + randn(2, 1).*[.1; .1]*dt;
%		y1(t, :) = x1';
%	end
%	[ms1, covs1, sigs1] = trajgpPred(Hyp, S.D(:, 1:3), S.D(:, 4:5), x0', diag([1e-6, 1e-6]), (10/dt/gpstep)+1, control_input, iscov, iicov, K1(:)', limit_state);
%
%	% Prediction with other starting state and same gain
%	x02 = (rand(1, 2)-.5)*2;
%	tm = [0:dt:10];
%	diffx_loc = @(t, x) diffx_int2(t, x, K);
%	%[t,y] = ode45(diffx_loc, tm, x0);
%	x2 = x02';
%	y2 = zeros(length(tm), 2);
%	y2(1, :) = x2';
%	for t = 2:length(tm)
%		x2 = x2 + diffx_loc(0, x2)*dt + randn(2, 1).*[.1; .1]*dt;
%		y2(t, :) = x2';
%	end
%	[ms2, covs2, sigs2] = trajgpPred(Hyp, S.D(:, 1:3), S.D(:, 4:5), x02', diag([1e-6, 1e-6]), (10/dt/gpstep)+1, control_input, iscov, iicov, K(:)', limit_state);

	figure(1)
	clf
	tm = [0:dt:10]';
	plot(tm, y, '--')
	title(['Same K, same x0, K=',num2str(K')])
	
	tm = [0:dtgp:10]';
	hold on
	plot(tm, ms(:, 1), 'b')
	plot(tm, ms(:, 1) + 2*sigs(:, 1), 'b-.')
	plot(tm, ms(:, 1) - 2*sigs(:, 1), 'b-.')
	plot(tm, ms(:, 2), 'r')
	plot(tm, ms(:, 2) + 2*sigs(:, 2), 'r-.')
	plot(tm, ms(:, 2) - 2*sigs(:, 2), 'r-.')

%	figure(4)
%	clf
%	tm = [0:dt:10]';
%	plot(tm, y1, '--')
%	title(['Altered K, K=',num2str(K1')])
%	tm = [0:dtgp:10]';
%	hold on
%	plot(tm, ms1(:, 1), 'b')
%	plot(tm, ms1(:, 1) + 2*sigs1(:, 1), 'b-.')
%	plot(tm, ms1(:, 1) - 2*sigs1(:, 1), 'b-.')
%	plot(tm, ms1(:, 2), 'r')
%	plot(tm, ms1(:, 2) + 2*sigs1(:, 2), 'r-.')
%	plot(tm, ms1(:, 2) - 2*sigs1(:, 2), 'r-.')
%
%	figure(5)
%	clf
%	tm = [0:dt:10]';
%	plot(tm, y2, '--')
%	title(['Altered x0, K=',num2str(K')])
%	tm = [0:dtgp:10]';
%	hold on
%	plot(tm, ms2(:, 1), 'b')
%	plot(tm, ms2(:, 1) + 2*sigs2(:, 1), 'b-.')
%	plot(tm, ms2(:, 1) - 2*sigs2(:, 1), 'b-.')
%	plot(tm, ms2(:, 2), 'r')
%	plot(tm, ms2(:, 2) + 2*sigs2(:, 2), 'r-.')
%	plot(tm, ms2(:, 2) - 2*sigs2(:, 2), 'r-.')

%	figure(2)
%	clf
%	sample_per_round = size(allData, 1)/e;
%	for i = 1:e
%		ix = [((i-1)*sample_per_round+1):(i*sample_per_round)];
%		plot3(allData(ix, 1), allData(ix, 2), allData(ix, 3), 'b'), hold on,
%	end
%	plot3(S.D(:, 1), S.D(:, 2), S.D(:, 3), 'k*'), grid on
%	xlabel('x1'), ylabel('x2'), zlabel('u')
%
%	figure(3)
%	plot(sort(ll))
%	title('likelihood values of training outputs under the current model')
%	xlabel(['individual loglikelihoods: ',num2str(logl)])

	pause

end
