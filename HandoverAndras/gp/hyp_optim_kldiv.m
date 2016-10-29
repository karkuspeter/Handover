function [f, df] = hyp_optim_kldiv(loghyp, x, y)

ix = randperm(length(y));

xtrain = x(ix(1:2:end), :); ytrain = y(ix(1:2:end), :);
xtest = x(ix(2:2:end), :); ytest = y(ix(2:2:end), :);


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
        
        % grad w
        for j = 1:d
            Wloc = zeros(d);
            Wloc(j, j) = w(j)^-3;
            gradK = k.* maha(xtrain, xtest(ii), Wloc);
            grad_mean(j, ii) = (gradK' - beta*gradLoc(:, :, j)) * alpha;
            grad_var(j, ii) = - 2*gradK' * beta' + beta*gradLoc(:, :, j)* beta';
        end
        
        % grad sigf
        grad_mean(d+1, ii) = 2/sigf*(k' - beta*K)*alpha;
        grad_var(d+1, ii) = 2*sigf - 2*2/sigf * k'*beta' + 2/sigf * beta * K * beta';
        
        
        % grad sign
        grad_mean(d+2, ii) = - 2*sign * beta * alpha;
        grad_var(d+2, ii) =   2*sign * (beta * beta') + 2*sign;
        
    end
    
    grad = grad_mean * (-2*pred_err./pred_var) + grad_var * ( -(pred_err./pred_var).^2 + 1./pred_var);
        
    f = f - mean( -.5*log(2*pi) - .5*log(pred_var) - .5 * pred_err.^2./pred_var);
    df(:, i) = .5*grad.*hyp(:, i)/nTest;
    
    disp([f, df(:, i)' hyp(:, i)']),
    
    
    [xsorted, ix] = sort(xtest);
%     figure(3), clf, plot(xsorted, pred_mean(ix), 'bo'), hold on, 
%     plot(xsorted, ytest(ix), 'r+'), drawnow
    
end

if imag(f) ~= 0 
 f = NaN;
end