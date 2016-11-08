close all, clear all

dt = .1;
tend = 10;
intstep = 10;
time = 0:dt:tend;

%matlabpool local 2
%matlabpool close
thres = .0001; N = 100; kernel_pars = [.5 , .5, .5, .5, .5, .2, .2, .2, .2]; 
S = {};

Qh = diag([1, 1, 1, 1, .1]);

logh = repmat(log([.5, .5, .5, .5, .05, 1., .3]), 4, 1);
iscov = @(S, K) K'*S;
iicov = @(S, K) K'*S*K;

w_pred = [.3, .3, .3, .3];
Wpred = diag(w_pred.^-2);

Vhyp = log([1, 1, 1, 1, 1, .1]);

Res = [];

K = randn(4, 1)*.1;

for e = 1:10
	xsave = [];
	xobssave = [];
	x = [.0, .0, .0, .1];
	sigs = [.1, .1, .1, .1];
	xobs = x + randn(1, 4).*sigs*dt;
	xsave = x;

	go = 1;
	counter = 1;

	while go

		if e > 1
			v = exp(-.5*diag( (Res(:, 1:4) - repmat(xobs(:)', size(Res, 1), 1))*Wpred*(Res(:, 1:4) - repmat(xobs(:)', size(Res, 1), 1))'));
			q = Res(:, 14).*v(:);
			K = Res(:, 5:8)'*q(:)/sum(q);
		end

		u = xobs*K;
		cartpole_w_input = @(t, y) cartpole(t, y, u);

		tm = [0:dt/intstep:dt];
		[t, xint] = ode45(cartpole_w_input, tm, x);
		xint(:, 4) = mod(xint(:, 4)+pi, 2*pi)-pi;
		x = xint(end, :);
		xobs = x+randn(1, 4).*sigs*dt;
		xobssave = [xobssave; [xobs, u]];
		xsave = [xsave; xint(2:end, :)];

		if counter > 1
			dnew = [xobssave(end-1, :), xobs - xobssave(end-1, 1:4)]; 
			[S, update, ix_update] = incr_online_sparsification(dnew, thres, N, kernel_pars, S, 0);
		end

		counter = counter + 1;

		if (counter > 500) || (abs(xobs(1)) > 1)
			go = 0;
		end

	end

	plot(xsave)
	pause

	for i=1:4
		options = optimset('GradObj', 'on');
		[hyp_opt, fval] = fminunc(@(hyp) hyper_optim(hyp, S.D(:, 1:5), S.D(:, 5+i), 0), logh(i, :), options);
		logh(i, :) = hyp_opt;
		disp(exp(hyp_opt))
	end

	input = S.D(:, 1:5);
	target = S.D(:, 6:end);


	for p=1:size(input, 1)

		for j = 1:1
			K = randn(4, 1)*5;
			m = input(p, :);
			e(1:4) = m(1:4) + randn(1, 4)*.1;
			m(end) =m(1:4)*K;
			m = m';
			s = .000001*ones(5, 5);

			[msave, ssave] = trajgpPred(logh, input, target, m, s, 50, iscov, iicov, K);

			l = 0;
			for kk = 1:size(ssave, 3)
				%l = l + trace(Qh*[ssave(:, :, kk), ssave(:, :, kk)*K; K'*ssave(:, :, kk), K'*ssave(:, :, kk)*K]) + [msave(kk, :), msave(kk, :)*K]*Qh*[msave(kk, :), msave(kk, :)*K]';
			end
			l = l/size(ssave, 3);

			Res = [Res; [m(1:4)', K', l, NaN, NaN, NaN, NaN, NaN]];

		end
	
	end

	% Hyperparameter optimization for value function prediction
	Vhypp = Vhyp;
	try
		options = optimset('GradObj', 'on');
		[Vhyp, fval] = fminunc(@(hyp) hyper_optim(hyp, Res(:, 1:4), Res(:, 9), 0), Vhyp, options);
		exp(Vhyp)
	catch
		disp('Vhyper optim failed')
		Vhyp = Vhypp;
	end


	% Predict the value function for each state
		Whyp = diag(exp(Vhyp(1:4)).^-2);
		VK = exp(Vhyp(5))^2*exp(-.5*maha(Res(:, 1:4), Res(:, 1:4), Whyp));
		ViKy = inv( VK + eye(size(Res, 1))*exp(Vhyp(6))^2);
		
		for i = 1:size(Res, 1)
			k = VK(i, :);
			Res(i, 10) = k*ViKy*Res(:, 9); % mean value func
			Res(i, 11) = exp(Vhyp(5))^2 - k*ViKy*k'; % var value func
			Res(i, 12) = Res(i, 10) - Res(i, 9); % advantage func (V - l)
			Res(i, 13) = -log(1/99)/.5/Res(i, 10); % eta for R = 1/(1+exp(-eta*A))
			Res(i, 14) = 1/(1+exp(-Res(i, 13)*Res(i, 12))); % transformed rewards

		end

	Res = Res(Res(:, 14) > 1e-5, :)

end
