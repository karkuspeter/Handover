clear all, close all, try reset(gpuDevice); end

addpath('/home/kupcsik/Dropbox/workonline/progs/matlab/gpml'); startup

load trajs;
input = [];
target = [];
dummy1 = [];
dummy2 = [];
for i = 1:10
    input = [input; trajs(1:end-1, (i-1)*12 + (1:12))];
    target = [target; diff( trajs(:, (i-1)*12 + [2 4 6 8]))];
    
    
    dummy2 = [dummy2; trajs(2:end, (i-1)*12 + (1:12))];
end

dummy2 = dummy2(:, 1:8);
dummy2 = dummy2 + bsxfun(@times, randn(size(dummy2, 1), 8), repmat([.001, .02 ]/5, 1, 4));

% transforming stuff
AB = input\(dummy2(:, 1:8)-input(:, 1:8));
target = dummy2(:, 1:8)-input(:, 1:8) - input*AB;


[n, d] = size(input);
e = size(target, 2);
hypInitX = std(input)/5;
logHypInit(1:d, :) = log(repmat(hypInitX(:), 1, e));
for i = 1:e
    logHypInit(d+1, i) = log(std(target(:, i)));
    logHypInit(d+2, i) = log(exp(logHypInit(d+1, i))/10);
end



i =3;
switch i 
    case 1 % gradient check
        
        j = 1;
        hyp_init = logHypInit(:, j);
        covfunc = @covSEard;
         likfunc = @likGauss;
         hyp.cov = hyp_init(1:end-1);
         hyp.lik = hyp_init(end);
         hyp.mean = []; 
         
         %first llhood and gradient (gpml)
         [f, df] = gp(hyp, @infExact, [], covfunc, likfunc, input, target(:, j));
         disp(['GPML f: ', num2str(f), ', df: ', num2str([df.cov(:)', df.lik])])          
         
         % first llhood and gradient (own)
         [f, g] = hyper_optim_GPUoptim_doubleOnly(hyp_init, input, target(:, j), 0);
         disp(['GPUd f: ', num2str(f), ', df: ', num2str(g(:)')])
         
         % first llhood and gradient (own)
         [f, g] = hyper_optim_GPUoptim(hyp_init, input, target(:, j), 0);
         disp(['GPUs f: ', num2str(f), ', df: ', num2str(g(:)')])
         
         
         
    case 2 % optimized
        for j = 1:4
            tic; 
            [hyp_opt, fval, lines] = minimizeGPU(logHypInit(:, j), @(hyp) hyper_optim_CPUoptim2(hyp, input, target(:, j), 0), -200);
            finish = toc;
            disp(['Optimized CPU version: ', num2str(finish), ' sec, fevals: ', num2str(lines), ', fevals/second: ', num2str(lines/finish)]); 
        end
    case 3 % gpu optimized
        profile on
        hyp_opt = getFullGPModelGPU(input, target, -200, 0);
%         hyp_opt = getFullGPModelGPUpso(input, target, 30, -100);
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
