function [ms, covs, sigs] = trajGpPredPseudo(logh, input, pinput, target, m, s, steps, control_input, iscov, iicov, K)
% Trajectory prediction with a GP model
% (Note: current assumption is that input state/output state dimensionality is equal)
%
% Inputs:
%	logh: the log of hyperparameters [D+2, E], (E: target/state dimensionality
% 												D: input+state dimensionality)
%   input: the training input [n, D] (n: number of training samples)
%	pinput: the pseudo  input [np, D]
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
%		output: column vector of control input [control_dim, 1]
% 	iscov: input-state covaraince [func] (inputs: {state cov, control param})
%		outputs: [control_dim, E] dimensional matrix
% 	iicov: input covaraince [func] (inputs: {state cov, control param})
% 	K: controller parameters [control_dim, p] or its transposed, depends on control_input, iscov, iicov
%	
% Outputs:
%	ms: the state mean predictions [steps, E+control_dim]
%	covs: the state covariance predictions [E+control_dim, E+control_dim, steps]
%	sigs: the state standard deviation [steps, E+control_dim]

if size(steps, 1) > 1
	tracking = 1;
	ref = steps;
	steps = size(ref, 1);
else
	tracking = 0;
end

[n, D] = size(input);
[n, E] = size(target);

control_dim = D-E;

logh_col = reshape(logh', E*(D+2), 1);

% Init mean and covariance
m = [m; zeros(control_dim, 1)];
sd = zeros(D, D);
sd(1:E, 1:E) = s;
s = sd;

% Update the first input distribution
if tracking
	u = control_input(m(1:E, 1), K, ref(1, :));
else
	u = control_input(m(1:E, 1), K);
end

m(E+1:end, 1) = u;
s(E+1:end, 1:E) = iscov(s(1:E, 1:E), K);
s(1:E, E+1:end) = iscov(s(1:E, 1:E), K)';
s(E+1:end, E+1:end) = iicov(s(1:E, 1:E), K);

ms(1, :) = m';
covs(:, :, 1) = s;
sigs(1, :) = diag(s)'.^.5;

% Some parameters for the trajectory prediction
gpmodel.hyp = logh;
gpmodel.inputs = input;
gpmodel.target = target;
gpmodel.induce = pinput;
shit = diag(exp(2*logh(end,:)));
for i=2:steps

	[M, S, V] = gp1(gpmodel, m, s);
	V = s*V;
    S = S+shit; % noisy measurement
	m(1:E, 1) = m(1:E, 1) + M; 

	if ~tracking
		m(E+1:end, 1) = control_input(m(1:E, 1), K);
	else
		m(E+1:end, 1) = control_input(m(1:E, 1), K, ref(i, :));
	end
	dummy = s(1:E,1:E) + S + V(1:E, :) + V(1:E, :)';
	s = [dummy, iscov(dummy, K)'; iscov(dummy, K), iicov(dummy, K)];

	

	ms(i, :) = m';
	covs(:, :, i) = s;
	sigs(i, :) = diag(s)'.^.5;

end
