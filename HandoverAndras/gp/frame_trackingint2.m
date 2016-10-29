clear all, close all

addpath('../gpml')
addpath('../SPGP_dist')
startup % of GPML

% Prediction parameters
control_input = @(x, k, r) (x(:)'-r(:)')*k(:);
iscov = @(S, K) K(:)'*S;
iicov = @(S, K) K(:)'*S*K(:);
limit_state = @(x) x;

sig_noise = [.1, .2];
K = -[4, 2];
x0 = [0, 0];
dt = .01;
tm = [0:dt:10];

diffx_int2 = @(t, x, r, k) [x(2), (x(:)'-r)*k(:)]';
diffx_locr = @(t, x, r) diffx_int2(t, x, r, K);
ref = [0, 1];
diffx_loc = @(t, x) diffx_locr(t, x, ref);

% mean trajectory
tm = [0:dt:10];
x = x0';
xsave(1, :) = x';
yr = zeros(length(tm), 2);
var = zeros(2, 2, length(tm));
sig = zeros(length(tm), 2);
var(:, :, 1) = diag(sig_noise.^2*dt);
sig(1, :) = sig_noise*dt^.5;
F = [0 1 0; 0 0 1]*dt;
for t = 2:length(tm)
	x = x + diffx_loc(0, x)*dt;
	yr(t, :) = x';
	var(:, :, t) = diag((sig_noise.^2)*dt^2*t) + F*[var(:, :, t-1), var(:, :, t-1)*K'; K*var(:, :, t-1), K*var(:, :, t-1)*K']*F';
	sig(t, :) = diag(var(:, :, t)).^.5;
end

tm = [0:dt:10];
xsave = [];
xsaveex = [];
usave = [];
notrajs = 10;
for e = 1:notrajs
	x = x0';
	y = zeros(length(tm), 2);
	for t = 2:length(tm)
		x = x + diffx_loc(0, x)*dt + randn(2, 1).*sig_noise'*dt*(t*dt)^.5;
		y(t, :) = x';
	end
	xsave = [xsave; y];
	xsaveex = [xsaveex y];
	usave = [usave; (y-repmat(ref, size(y,1), 1))*K'];
end



trainx = [xsave, usave];
train_y = [];
train_x = [];
for i = 1:notrajs
	train_y = [train_y; diff(trainx((i-1)*length(tm)+2:i*length(tm), 1:2))];
	train_x = [train_x; trainx( (i-1)*length(tm)+2:i*length(tm)-1, :)];
end

[hyp, px] = hyper_pseudo_commonInputs(train_x, train_y, 50, 200); 

[ms, covs, sigs] = trajGpPredPseudo(hyp, train_x, px, train_y, x0', diag([1e-6, 1e-6]), repmat(ref, length(tm), 1), control_input, iscov, iicov, K(:)', limit_state);

figure, 

plot(tm, ms(:, 1), 'r'), hold on,
plot(tm, ms(:, 1)+2*sigs(:, 1), 'r--')
plot(tm, ms(:, 1)-2*sigs(:, 1), 'r--')
plot(tm, ms(:, 2), 'b')
plot(tm, ms(:, 2)+2*sigs(:, 2), 'b--')
plot(tm, ms(:, 2)-2*sigs(:, 2), 'b--')

plot(tm, yr(:, 1), 'r', 'LineWidth', 2), hold on,
plot(tm, yr(:, 1) + 2*sig(:, 1), 'r--', 'LineWidth', 2)
plot(tm, yr(:, 1) - 2*sig(:, 1), 'r--', 'LineWidth', 2)
plot(tm, yr(:, 2), 'b', 'LineWidth', 2), hold on,
plot(tm, yr(:, 2) + 2*sig(:, 2), 'b--', 'LineWidth', 2)
plot(tm, yr(:, 2) - 2*sig(:, 2), 'b--', 'LineWidth', 2)

ylabel('GP/Real(thick)')

[kldiv, ikldiv] = KLdivTraj(ms, covs, yr, var);
figure
semilogy(tm, kldiv, 'r'), hold on
semilogy(tm, ikldiv, 'b')
legend('GP||Real','Real||GP')
grid

