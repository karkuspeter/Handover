clear all, close all

N = 200;
x = 10*(rand(N, 1)-.5);
x = sort(x);
y = sinc(x) + randn(N, 1)*.1;

S = {};
for i = 1:length(x)
	dnew = [x(i) y(i)]; thres = .01; N = 15; kernel_pars = [.4, .2]; 
	[S, update, ix_update] = incr_online_sparsification(dnew, thres, N, kernel_pars, S, 0);
end

hyp = log([.4, 1, .1]);

options = optimset('GradObj', 'on');
[hyp_opt, fval] = fminunc(@(hyp) hyper_optim(hyp, S.D(:, 1), S.D(:, 2), 0), hyp, options)
disp(exp(hyp_opt))

xsamp = S.D(:, 1);
ysamp = S.D(:, 2);

w = exp(hyp_opt(1));
sigf = exp(hyp_opt(2));
signs = exp(hyp_opt(3));


K = sigf^2*exp(-.5*maha(xsamp, xsamp, 1/w^2));
iKy = inv(K + eye(length(xsamp))*signs^2);

predx = zeros(length(x), 1);
varx = zeros(length(x), 1);
x = sort(x);
for i = 1:length(x)
	ker = sigf^2*exp(-.5/w^2*(xsamp-x(i)).^2);
	predx(i) = ker(:)'*iKy*ysamp(:);
	varx(i) = sigf^2 - ker(:)'*iKy*ker(:);
end

plot(x, predx, 'k', 'Linewidth', 2); hold on
plot(x, predx + varx.^.5, 'r');
plot(x, predx - varx.^.5, 'r');


plot(sort(x), sinc(sort(x))), 
plot(S.D(:, 1), S.D(:, 2), 'k*')
