function e = evidence_plgp(loghyp, x, prefs, ridge)

sig = exp(loghyp(1));
w = exp(loghyp(2:end));

K = exp(-.5 * maha(x(:), x(:), w));
K = K + eye(size(K))*ridge;

f = randn(length(x), 1);
optfun = @(ff) grad_fmap(ff, prefs, K, sig);

% options_loc = optimoptions('fminunc','GradObj','on', 'Hessian', 'on');

fmap = nr_plgp(f, prefs, K, sig);

[S, dS, ddS, beta] = optfun(fmap);

GammaMap = ddS - eye(size(K))/(K + eye(size(K))*ridge);

e = S  + .5* log(det(K*ddS + eye(size(K)) ));