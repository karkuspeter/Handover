clear all, close all

x = -15:.011:15;
y = sin(x)./x + randn(1, length(x))*.1;
yy = y;

% plot(x, y, '*')


m = 15;
[hyp, hyp_orig px] = hyper_pseudo([x', sin(x)'],  y', m, 100);

[mu1,s1,t_train,t_test] = spgp_pred(y',[x', sin(x)'], px, [x', sin(x)'],hyp_orig);

figure, plot(x, mu1), hold on, plot(x, mu1+2*s1.^.5, 'b--'), plot(x, mu1-2*s1.^.5, 'b--')
plot(x, sin(x)./x, 'r')

% 
% 
% profile on
% % Precomputation
% w = exp(hyp(1:end-2));
% W = diag(w.^-2);
% sign_sq = exp(hyp(end))^2;
% sigf_sq = exp(hyp(end-1))^2;
% n = length(x);
% tic
% K = sigf_sq*exp(-.5*maha(px, px, W));
% L = chol(K)';
% x1 = x'.*repmat((1./w)',n,1);
% x2 = px.*repmat((1./w)',m,1);
% K = -2*x1*x2' + bsxfun(@plus, sum(x2.*x2,2)', sum(x1.*x1,2));
% clear x1 x2
% K = sigf_sq*exp(-.5*K');
% V = L\K;
% ep = 1 + (sigf_sq*ones(n, 1)-(ones(1, m)*V.^2)')/sign_sq;
% V = V./repmat(sqrt(ep)',m,1); 
% y = y'./sqrt(ep);
% Lm = chol(sign_sq*eye(m) + V*V')';
% bet = Lm\(V*y);
% clear V
% t_train_mine = toc;
% 
% [L2, Lm2, bet2] = getLLmbet(x', yy', px, hyp);
% % Prediction
% tic
% x1 = x'.*repmat((1./w)',n,1);
% x2 = px.*repmat((1./w)',m,1);
% K = -2*x1*x2' + bsxfun(@plus, sum(x2.*x2,2)', sum(x1.*x1,2));
% K = sigf_sq*exp(-.5*K');
% lst = L\K;
% clear K
% lmst = Lm\lst;
% mu2 = (bet'*lmst)';
% s2 = sigf_sq*ones(n, 1) - (ones(1, m)*lst.^2)' + (sign_sq*ones(1, m)*lmst.^2)';
% t_pred_mine = toc;
% 
% figure, 
% plot(mu1, 'r'), hold on, 
% plot(mu1+2*s1.^.5, 'r--')
% plot(mu1-2*s1.^.5, 'r--')
% plot(mu2, 'b'), hold on, 
% plot(mu2+2*s2.^.5, 'b--')
% plot(mu2-2*s2.^.5, 'b--')
% profile viewer