clear all, close all

rew = @(x, y) exp(-.5*(x.^2/.4^2 + y.^2/.4^2));
rew = @(x, y) -(20 + x.^2 + y.^2 - 10*(cos(2*pi*x) + cos(2*pi*y)));

[X, Y] = meshgrid(-1:.01:1, -1:.01:1);
Z = rew(X, Y);
figure()
contourf(X, Y, Z, 15), colorbar, hold on,
%rew = @(x, y) input(['Reward for [', num2str(x), ',', num2str(y),']?: ']);

N = 30;
episodes = 50;
epsilon = .5;
eta = .5;

mn = [0; -.5];
cn = diag([.5, .5].^2);

savem = [mn', diag(cn)'.^.5, cn(2, 1)/prod(diag(cn).^.5), eta, 0, 0, 0, 0];

options = optimset('Algorithm','active-set');
options = optimset(options, 'GradObj','on');
options = optimset(options, 'Display', 'off');
S = [];
r = [];

for e = 1:episodes

	Snew = repmat(mn', N, 1) + randn(N, 2).*repmat(diag(cn)'.^.5, N, 1);
	%q = exp(-.5*diag((Snew - repmat(mn', N, 1))*inv(cn)*(Snew - repmat(mn', N, 1))'));
	%q = q/sum(q);
	q = ones(N, 1)/N;
    
        rnew = rew(Snew(:, 1), Snew(:, 2));
    
	%S = [S; Snew];
	%r = [r; rnew];
	
	S = Snew;
	r = rnew;

	% With Gradient
	objfun = @(eta) reps_dual(eta, r, q, epsilon);
	logeta = fmincon(objfun, log(abs(mean(r))), -1, -.01, [], [], [], [], [], options);
    
    eta = exp(logeta);

	p = exp(r/eta)/sum(exp(r/eta));
	Dkl = sum(p.*log(p./q));
	iDkl = sum(q.*log(q./p));
	mn = S'*p;
	cn = (S-repmat(mn', size(S, 1), 1))'*((S-repmat(mn', size(S, 1), 1)).*repmat(p, 1, 2));

	rexp = r(:)'*p;
	rvar = (r(:)-rexp)'*((r(:)-rexp).*p);

	savem = [savem; mn', diag(cn)'.^.5, cn(2, 1)/prod(diag(cn).^.5), eta, rexp, Dkl, iDkl, rvar];
end

for i = 1:(size(savem, 1)-1)
	th = .5*atan(2*savem(i+1, 5) * savem(i+1, 3) * savem(i+1, 4)/(savem(i+1,3)^2 - savem(i+1, 4)^2));
	[ex, ey] = calculateEllipse(savem(i+1, 1), savem(i+1, 2), savem(i+1, 3), savem(i+1, 4),th *180/pi); 
	plot([savem(i, 1); savem(i+1, 1)], [savem(i, 2), savem(i+1, 2)], 'w', 'Linewidth', 2), hold on;
	plot(ex, ey, 'w--')
end

figure()
subplot(4,2,1), plot(savem(1:end, 1:2)), ylabel('mean')
subplot(4,2,2), plot(savem(2:end, 7)), ylabel('E[R]')
subplot(4,2,3), plot(savem(1:end, 3:4)), ylabel('std')
subplot(4,2,4), plot(savem(2:end, 6)), ylabel('eta')
subplot(4,2,5), plot(savem(2:end, 8)), ylabel('Dkl')
subplot(4,2,6), plot(savem(1:end, 9)), ylabel('iDkl')
subplot(4,2,7), plot(savem(2:end, 10).^.5), ylabel('Var[R]^{0.5}')
