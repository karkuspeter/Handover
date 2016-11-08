function [g, dg] = reps_dual_state(etatheta, epsilon, R, Phi)

eta = exp(etatheta(1));
theta = etatheta(2:end);

N = length(R);
V = Phi*theta(:);

MR = max(R(:)-V);
Zhat = 1/N* exp((R(:)-V - MR)/eta);

g = eta*log(sum(Zhat)) + MR + eta*epsilon + mean(Phi, 1)*theta(:);
dgeta = epsilon + log(sum(Zhat)) + MR/eta - sum(Zhat.*(R(:)-V))/eta/sum(Zhat);
dgeta= dgeta *eta;
dgtheta = mean(Phi, 1)' - Phi'*Zhat/sum(Zhat);

dg = [dgeta, dgtheta'];

% [g, dg]

% if any(isnan([g, dg]))
%     warning('nan g or dg in contextual reps dual')
% 	
% end

