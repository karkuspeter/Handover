function [g, dg] = reps_dual(eta, r, epsilon)

a = r / eta - log(length(r));
A = max(a);
logSumExp = A + log(sum(exp(a - A)));

g = eta * (epsilon + logSumExp);
dg = epsilon + logSumExp - sum(exp(r/eta).*r) / length(r) / exp(logSumExp) / eta;