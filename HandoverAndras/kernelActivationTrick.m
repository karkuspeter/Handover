function w_opt = kernelActivationTrick(samples, activations)

w = medianTrick(samples, 1);
w_opt = zeros(length(activations), length(w));
options = optimoptions('fminunc', 'Algorithm','trust-region','GradObj','on','Hessian', 'off', 'MaxFunEvals', 100, 'TolFun', 1e-100, 'Display', 'off');
for i = 1:length(activations)
    
    optfun = @(logwScaler) kernelActivationOptFun(logwScaler, activations(i), samples, w);
    logScaler_opt = fminunc(optfun, log(1), options);
    
    w_opt(i, :) = exp(logScaler_opt) * w(:)';
end
    

    
    