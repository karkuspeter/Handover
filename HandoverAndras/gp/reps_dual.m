function [g, dg] = reps_dual(eta, r, q, epsilon);

Z = exp(r/eta).*q;
g = eta * log(sum(Z)) + eta*epsilon;
dg = epsilon + log(sum(Z)) - sum(Z.*r)/(eta * sum(Z));
