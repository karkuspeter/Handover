function [f, df] = hyp_optim_kldiv_numerical(loghyp, x, y)


func = @(par) getLikelihood(par, x, y);
df = numericalGradient(func, loghyp);

f = getLikelihood(loghyp, x, y);

disp([f, df(:)' exp(loghyp(:))'])

end

function f = getLikelihood(loghyp, x, y)
xtrain = x((1:2:end), :); ytrain = y((1:2:end), :);
xtest = x((2:2:end), :); ytest = y((2:2:end), :);


e = size(y, 2);
[nTrain, d] = size(xtrain);
nTest = size(xtest, 1);

if and(or(size(loghyp, 2) == 1, size(loghyp, 1) == 1), e > 1)
    loghyp = reshape(loghyp(:), d+2, e);
end

hyp = exp(loghyp);

df = zeros(d+2, e);
f = 0;

for i = 1:e
    sigf = hyp(end-1,i); sign = hyp(end, i); w = hyp(1:end-2, i);
    K = sigf^2 * exp(-.5*maha(xtrain, xtrain, diag(w.^-2)));
    L = chol(K+eye(nTrain)*sign^2);
    
    alpha = eye(nTrain)/L/L'*ytrain(:, i);
    
    grad_mean = zeros(d+2, nTest);
    grad_var = zeros(d+2, nTest);
    
    for j = 1:d
        Wloc = zeros(d);
        Wloc(j, j) = w(j)^-3;
        gradLoc(:, :, j) = K.*maha(xtrain, xtrain, Wloc);
    end
    
    pred_mean = zeros(nTest, 1);
    pred_var = zeros(nTest, 1);
    pred_err = zeros(nTest, 1);
    
    for ii = 1:nTest
        
        k = sigf^2 * exp(-.5*maha(xtrain, xtest(ii), diag(w.^-2)));
        beta = k'/L/L';
        
        pred_mean(ii) = k' * alpha;
        pred_err(ii) = ytest(ii, i) - pred_mean(ii);
        pred_var(ii) = sigf^2 - beta*k + sign^2;
    end
end
dummy =  -.5*log(2*pi) - .5*log(pred_var) - .5 * pred_err.^2./pred_var;

% figure(2), clf, plot(sort(-dummy));
f = f - sum(dummy);

% f = log(sum(pred_err.^2.*pred_var));
end

