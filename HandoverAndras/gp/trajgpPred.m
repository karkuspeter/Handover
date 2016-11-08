function [ms, covs, sigs] = trajgpPred(logh, input, target, m, s, steps, control_input, iscov, iicov, K, limit_state)
% Trajectory prediction with a GP model
%
% Inputs:
%	logh: the log of hyperparameters [E, D+2], (E: target/state dimensionality
% 												D: input+state dimensionality)
%   input: the training input [n, D] (n: number of training samples)
%	target: the training traget [n, E]
% 	m: initial mean [E, 1]
% 	s: initial covariance [E, E]
% 	steps: either
%        number of timesteps to predict [1]
%				--- OR ---
%		 reference of states for a tracking problem [steps, E]
%  	control_input: mapping from states to control input [func] 
%	  	inputs: either
%  			{state, control param} 
%				---- OR ---
%		    {state, control param, reference}
% 	iscov: input-state covaraince [func] (inputs: {state cov, control param})
% 	iicov: input covaraince [func] (inputs: {state cov, control param})
% 	K: controller parameters [1, p]
%   limit_state: manipulates states after mean update [func] (input: {state})
%	
% Outputs:
%	ms: the state mean predictions [steps, E]
%	covs: the state covariance predictions [E, E, steps]
%	sigs: the state standard deviation [steps, E]

if size(steps, 1) > 1
	trackig = 1;
	ref = steps;
	steps = size(ref, 1);
else
	tracking = 0;
end

[n, D] = size(input);
[n, E] = size(target);

control_dim = D-E;

msave = zeros(steps, E);
ssave = zeros(E, E, steps);

logh_col = reshape(logh', E*(D+2), 1);

m = [m; zeros(control_dim, 1)];
sd = zeros(D, D);
sd(1:E, 1:E) = s;
s = sd;

for i=1:steps

	[M, S, V] = gpP0(logh_col, input, target, m, s);
	V = s*V;
    S = S+diag(exp(2*logh_col(end,:))); % noisy measurement
	
	m(1:E, 1) = m(1:E, 1) + M; 
	m(1:E, 1) = limit_state(m(1:E, 1)); % is it necessary?
	if ~tracking
		m(E+1:end, 1) = control_input(m(1:E, 1), K);
	else
		m(E+1:end, 1) = control_input(m(1:E, 1), K, ref(i, :));
	end
	dummy = S + V(1:E, :) + V(1:E, :)';
	s = s + [dummy, iscov(dummy, K)'; iscov(dummy, K), iicov(dummy, K)];

	ms(i, :) = m(1:E, 1)';
	covs(:, :, i) = s(1:E, 1:E);
	sigs(i, :) = diag(s(1:E, 1:E))'.^.5;

end
