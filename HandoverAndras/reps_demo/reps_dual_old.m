function [g, dg] = reps_dual(eta, r, q, epsilon)

Z = exp(r/eta).*q ;

a = r/eta - log(length(r));
A = max(a);

logSumExp = A + log(sum(exp(a - A))); 

g = eta * logSumExp + eta*epsilon;
% dg = epsilon + logSumExp - sum(Z.*r)/(eta * sum(Z));


% dg = epsilon + logSumExp - sum(Z.*r) * exp( - log(eta) - A - log(sum(exp(a - A))));
dg = epsilon + logSumExp - sum(Z.*r) * exp( - log(eta) - logSumExp);

if isinf(g)
    keyboard
end

if imag(logSumExp) ~= 0
    keyboard
end

% [g, dg]

