function [g, dg] = reps_dual_state(etatheta, epsilon, R, Phi)

eta = etatheta(1);
theta = etatheta(2:end);

N = length(R);
V = Phi*theta(:);
MR = max(R-V);
Zhat = 1/N* exp((R-V - MR)/eta);

g = eta*log(sum(Zhat)) + MR + eta*epsilon + mean(Phi, 1)*theta(:);
dgeta = epsilon + log(sum(Zhat)) + MR/eta - sum(Zhat.*(R-V))/eta/sum(Zhat);
dgtheta = mean(Phi, 1)' - Phi'*Zhat/sum(Zhat);

dg = [dgeta, dgtheta'];


%if any(isnan([g, dg]))
%	keyboard
%end

