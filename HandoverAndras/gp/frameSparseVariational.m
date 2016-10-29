clear all, close all

addpath('../gpml')
startup

%% Init
xorig = [-15:.1:-.1, .1:.1:15]';
yorig = sin(xorig)./xorig + randn(length(xorig), 1)*.1;
n = length(xorig);
e = 1;
d = 1;
m = 10;

hyp = zeros(d+2, e);
px = zeros(m, d);
hyp_init = [];
randix = randperm(n);
maxm = length(randix);
m = min(m, maxm);
xb_init = xorig(randix(1:m), :);

for i = 1:e
		hyp_init_marc = [hyp_init; log(std(xorig)'); % log lengthscales
			log(std(yorig(:, i))); % log sf
			log((std(yorig(:, i)))/25)]; % log sn
% 		hyp_init_mine = [hyp_init; log(1./((std(x)/2)'.^2)); % log 1/(lengthscales)^2
% 			2*log(std(y(:, i))); % log sf^2
% 			2*log((std(y(:, i)))/25)]; % log sn^2
end
w_init = hyp_init_marc;
    
%% Optimize
useGPU = 0;
useSinglePrecision = 0;
if useGPU
	reset(gpuDevice);
end
ix = randperm(length(xorig));
hypFull = getFullGPModel(xorig, yorig);

px = xorig(ix(1:m));
 
w_init = [w_init; px];
 
% w = minimize(w_init, 'Var_lik', -50, xorig, yorig, m, useGPU, useSinglePrecision);
% hyp = w(1:3);
% px = w(4:end);
% ixPseudo = ones(length(x), 1) == 0;
% ixPseudo(ix(1:m)) = ones(m, 1) == 1;
% wm = hyp_init_mine;
% wM = hyp_init_marc;
% wMpseudo = [hyp_init_marc; px];
% wDiff = [wFull - wm];
% for i = 1:1
% 
% 	% M-step
% %	tic
% %	wm = minimize(wm, 'Var_lik_hypOnly', -50, px, xorig, yorig, ixPseudo, useGPU, useSinglePrecision);
% %	disp(['hyp only mine: ', num2str(toc)])
% %	tic
% %	wMpseudo = minimize(wMpseudo, 'Var_lik_hyp', -50, xorig, px, yorig, useGPU, useSinglePrecision);
% %	disp(['hyp Marc: ', num2str(toc)])
% 	tic
% 	wMpseudo = minimize(wMpseudo, 'Var_lik', -50, xorig, yorig, m, useGPU, useSinglePrecision);
% 	disp(['hyp + pseudo Marc: ', num2str(toc)])
% 	% E-step
% 	%Fv = varSparseGPEstep(w, px, xorig, yorig, ixPseudo);
% 	%ixNotPseudo = find(ixPseudo == 0);
% 	%[FvMax, ixMax] = max(Fv);
% 
% 	%px = [px; xorig(ixNotPseudo(ixMax))];
% 	%py = [py; yorig(ixNotPseudo(ixMax), :)];
% 	%ixPseudo(ixNotPseudo(ixMax)) = 1==1;
% 
% %	wDiff = [wDiff, wFull - exp(w)];
% end
% 
% % extend pseudo input set to 40
% %ixNotPseudo = find(ixPseudo == 0);
% %m = size(px, 1);
% %ixRand = randperm(length(ixNotPseudo));
% %if m < 100
% %	ixExtra = ixRand(1:(100-m));
% %	px = [px; xorig(ixExtra)];
% %	py = [py; yorig(ixExtra, :)];
% %end
[hyp, pxnew] = hyper_variationalPseudo_commonInputs(xorig, yorig, m, 100, [])
[mustar, Astar] = variationalSparsePrecomputations(hyp, pxnew, xorig, yorig);

Knm = exp(hyp(2))* exp(-.5*maha(xorig, pxnew, exp(hyp(1))));
mupred = Knm*mustar;
stdpred = (ones(n, 1)*exp(hyp(3)) + ones(n, 1)*exp(hyp(2)) - bsxfun(@dot, Knm', Astar*Knm')').^.5;



%% Plot
x = xorig;
y = yorig;
figure, %plot(x, y, '*')
plot(x, sin(x)./x, 'k')
hold on, plot(px, ones(length(px), 1)*(-.6), 'x') %plot(px, py(:, 1), 'r+')

hold on, plot(pxnew, ones(length(pxnew), 1)*(-.5), 'r+') %plot(px, py(:, 1), 'r+')
axis([-15, 15, -.8, 1.2])

[meanFull, sigFull] = predictWithFullGPModel(hypFull, xorig, yorig, xorig);
plot(x, meanFull, 'r')
plot(x, meanFull+2*sigFull, 'r--')
plot(x, meanFull-2*sigFull, 'r--')

plot(x, mupred, 'b');
plot(x, mupred+2*stdpred, 'b--');
plot(x, mupred-2*stdpred, 'b--');



