clear all, close all


episodes = 500;
batch = 1;
maxmembers = 100;
startlearn = 30;
epsilon = .7;


mn = [0, 0, 0, 0]';
cn = diag([2, 2].^2);

etatheta = [1, 1, 1];
savem = [mn', diag(cn)'.^.5, zeros(1, 7)];

options = optimset('Algorithm','active-set');
options = optimset(options, 'GradObj','on');
options = optimset(options, 'Display', 'off');
%options = optimset(options, 'MaxFunEval', 1000);
%options = optimset(options, 'MaxIter', 1000);
Aconst = -eye(3); b = [-1e-3, Inf, Inf];

Phi = []; A = []; S = []; R = [];

for e = 1:episodes

	Snew = [];
	Rnew = [];
	Anew = [];

	beta = [mn(1:2)'; mn(3:4)'];

	for i = 1:batch

		pos = rand(1)*3+1; % in range(3, 4)
		Snew = [Snew; pos];
		phi = [pos; pos^2];

		thv = mvnrnd(beta*[1, pos]', cn, 1);
		th = thv(1);
		v = thv(2);
		Anew = [Anew; thv];
		dist = distance_travelled(th, v, 1);
		r = -(pos-dist)^2;
		Rnew = [Rnew; r];
	end

	Phinew = [Snew, Snew.^2];
	A = [A; Anew];
	Phi = [Phi; Phinew];
	R = [R; Rnew];
	S = [S; Snew];

	if length(R) >= maxmembers
		dummy = length(R) - maxmembers + 1;
		R = R(dummy:end);
		A = A(dummy:end, :);
		Phi = Phi(dummy:end, :);
		S = S(dummy:end);
	end

	
	if length(R) >= startlearn
		epsilon_loc = length(R)/maxmembers * epsilon;
		objfun = @(etatheta) reps_dual_state(etatheta, epsilon_loc, R, Phi);
		etatheta = fmincon(objfun, etatheta, Aconst, b, [], [], [], [], [], options);

		eta = etatheta(1); theta = etatheta(2:end);
		V = Phi*theta(:);
		p = exp((R-V)/eta)/sum(exp((R-V)/eta));
		q = 1/length(p)*ones(length(p), 1);
		Dkl = sum(p.*log(p./q));
		iDkl = sum(q.*log(q./p));

		Sex = [ones(size(S, 1), 1), S];
		wSex = repmat(p, 1, 2).*Sex;
		mn(1:2) = (wSex'*Sex)\wSex'*A(:, 1);
		mn(3:4) = (wSex'*Sex)\wSex'*A(:, 2);
		beta = [mn(1:2)'; mn(3:4)'];
		cn = ((A' - (beta*Sex')).*repmat(p(:)', 2, 1))*(A' - (beta*Sex'))';

		rexp = R(:)'*p;
		rvar = (R(:)-rexp)'*((R(:)-rexp).*p);

		savem = [savem; mn', diag(cn)'.^.5, eta, theta, rexp, rvar, Dkl, iDkl];
		[mean(Phi, 1); p'*Phi];
	end


end

plot_reps_state(savem)

