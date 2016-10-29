clear all, close all, try reset(gpuDevice); end

addpath('~/Dropbox/workonline/progs/matlab/gpml'); 
startup

sinc = @(x) sin(x)./x;

width = 30;

x = -width/2:.11:width/2;
y = sinc(x);

N = 100;
xsamp = width*(rand(N, 1)-.5);
ysamp = sinc(xsamp) + randn(N, 1)*.1;

figure, plot(x, y, 'r', 'LineWidth', 2) 
hold on, plot(xsamp, ysamp, 'k*')

hyp_init = log([.5, .3, .1])';
options = optimset('GradObj', 'on');

% [f, df] = hyp_optim_kldiv_numerical(hyp_init, xsamp, ysamp)

hyp_opt = minimize(hyp_init, @hyp_optim_kldiv, 200, xsamp(:), ysamp(:));
hyp = exp(hyp_opt);

% hyp = [2.3 .2 .1];
% hyp = [0.3570    0.6844    0.1137];

sigf = hyp(end-1); 
sign = hyp(end); 
w = hyp(1:end-2);
K = sigf^2 * exp(-.5*maha(xsamp, xsamp, diag(w.^-2)));

alpha = (K + eye(N)*sign^2)\ysamp;

for i = 1:length(x)
    k = sigf^2 * exp(-.5*maha(xsamp, x(i), diag(w.^-2)));
    k = k(:);

    pred_mean(i) = k' * alpha;
    pred_var(i) = sigf^2 - k'/(K + eye(N)*sign^2)*k + sign^2;
end

hold on, plot(x, pred_mean, 'k', 'LineWidth', 2), hold on,
hold on, plot(x, pred_mean + 2*pred_var.^.5, 'k--', 'LineWidth', 1)
hold on, plot(x, pred_mean - 2*pred_var.^.5, 'k--', 'LineWidth', 1)