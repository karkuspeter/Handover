function [actError, dError] = kernelActivationOptFun(logwScaler, targetAct, samples, w)

Sigma = exp(-.5* maha(samples, samples, diag((exp(logwScaler)*w).^-2)));

actError = (median(mean(Sigma, 2)) - targetAct)^2;

if nargout > 1
    optfun = @(lws) kernelActivationOptFun(lws, targetAct, samples, w);
    dError = numericalGradient(optfun, logwScaler);
end