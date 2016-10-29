clear all, close all, 

addpath('~/Dropbox/workonline/progs/matlab/gpml'); startup

sinc = @(x) sin(x)./x;
x = -50:.01:50;
y = sinc(x);

N = 150;
xsamp = 50*2*(rand(N, 1)-.5);
ysamp = sinc(xsamp) + randn(N, 1)*4e-2;

% load jani;
% xsamp = jani.xsamp;
% ysamp = jani.ysamp;


plot(x, y); hold on,
plot(xsamp, ysamp, 'o')


% finding the optimal hyperparams3
% initial guess
hyp_init = log([.4, 1, .1]);

options = optimset('GradObj', 'on');
profile on

i =1;
switch i 
    case 0 % gradient check
        covfunc = @covSEard;
         likfunc = @likGauss;
         hyp.cov = hyp_init(1:2);
         hyp.lik = hyp_init(3);
         hyp.mean = [];
         
         %first llhood and gradient (gpml)
         [f, df] = gp(hyp, @infExact, [], covfunc, likfunc, xsamp(:), ysamp(:));
         disp(['GPML f: ', num2str(f), ', df: ', num2str([df.cov, df.lik])])          
         
         % first llhood and gradient (own)
%          [f, g] = hyper_optim_GPUoptim_doubleOnly(hyp_init, xsamp(:), ysamp(:), 0);
%          disp(['GPU (double precision) f: ', num2str(f), ', df: ', num2str(g(:)')])
         
         % first llhood and gradient (own)
%          [f, g] = hyper_optim_GPUoptim(hyp_init, xsamp(:), ysamp(:), 0);
%          disp(['GPU (single precision) f: ', num2str(f), ', df: ', num2str(g(:)')])
         
         
         
    case 1 % standard gpml
        tic
         covfunc = @covSEard;
         likfunc = @likGauss;
         hyp.cov = hyp_init(1:2);
         hyp.lik = hyp_init(3);
         hyp.mean = [];
         
         %first gradient
         [f, df] = gp(hyp, @infExact, [], covfunc, likfunc, xsamp(:), ysamp(:));
         disp(['GPML f: ', num2str(f), ', df: ', num2str([df.cov, df.lik])])          
         
        hyp_opt = minimize(hyp, @gp, -200, @infExact, [], covfunc, likfunc, xsamp(:), ysamp(:));
        disp(exp(hyp_opt.cov)); disp(exp(hyp_opt.lik))
        disp(['GPML CPU version: ', num2str(toc), ' sec']); 
    case 2 % optimized
        tic; [hyp_opt, fval, lines] = minimizeGPU(hyp, @(hyp) hyper_optim_CPUoptim(hyp, xsamp(:), ysamp(:), 0), -200);
	finish = toc;
        disp(['Optimized CPU version: ', num2str(finish), ' sec, fevals: ', num2str(lines), ', fevals/second: ', num2str(lines/finish)]); 
    case 3 % gpu optimized
        hyp_opt = getFullGPModelGPU(xsamp(:), ysamp(:), -200);
	
end
        
profile viewer
disp(exp(hyp_opt))
% hyp_opt = exp(hyp_opt);
disp(['SNR: ', num2str(exp(hyp_opt(end-1))/exp(hyp_opt(end)))]);
% 
% if i == 3
%     xsamp = gather(xsamp);
%     ysamp = gather(ysamp);
%     hyp_opt = gather(hyp_opt);
% end
% 
% w = exp(hyp_opt(1));
% sigf = exp(hyp_opt(2));
% signs = exp(hyp_opt(3));
% 
% 
% K = sigf^2*exp(-.5*maha(xsamp, xsamp, 1/w^2));
% iKy = inv(K + eye(length(xsamp))*signs^2);
% 
% predx = zeros(length(x), 1);
% varx = zeros(length(x), 1);
% for i = 1:length(x)
% 	ker = sigf^2*exp(-.5/w^2*(xsamp-x(i)).^2);
% 	predx(i) = ker(:)'*iKy*ysamp(:);
% 	varx(i) = sigf^2 - ker(:)'*iKy*ker(:);
% end
% 
% plot(x, predx, 'k', 'Linewidth', 2);
% plot(x, predx + varx.^.5, 'r');
% plot(x, predx - varx.^.5, 'r');
% ws = .01:.01:1;
% sigfs = .1:.02:2;
% 
% logl = zeros(length(ws), length(sigfs));
% xx = zeros(length(ws), length(sigfs));
% yy = zeros(length(ws), length(sigfs));
% 
% 
% for i=1:length(ws) % weights
% 	for j = 1:length(sigfs)
% 		loghyp = log([ws(i), sigfs(j), signs]);
% 		xx(i, j) = ws(i);
% 		yy(i, j) = sigfs(j);
% 		[f, g] = hyper_optim(loghyp, xsamp(:), ysamp(:), 0);
% 		logl(i, j) = -f;
% 	end
% end
% 
% figure();
% contour(xx, yy, logl, 20); hold on
% ylabel('sigf'); xlabel('w');
% title(['The log-likelihood with sign = ', num2str(signs)])
% colorbar
