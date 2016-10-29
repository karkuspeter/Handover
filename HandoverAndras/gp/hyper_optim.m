function [f, g] = hyper_optim(loghyp, input, target, verb)


% hyp: the hyperparameters w, sigf, sign
% input: n x d matrix
% target: n x 1 vector
% verb: verbosity

[n, d] = size(input);

w = exp(loghyp(1:d));
sigf = exp(loghyp(d+1));
signs = exp(loghyp(d+2));

W = diag(w.^-2);
K = sigf^2 * exp(-.5*maha(input, input, W));
Ky = K + eye(n)*signs^2;

iKy = inv(Ky);
dKy = det(Ky);

f = -(-.5*target'*iKy*target-.5*log(dKy) - n/2*log(2*pi)) ;
% for i = 1:d
% 	f = f + (w(i)/std(input(:, i))/100)^10;
% end
% maxSNR = 1000;
% f = f + (sigf/signs/maxSNR)^10;

alpha = iKy*target;
dummy = alpha*alpha'-iKy;

g = zeros(length(loghyp), 1);

% dK / dsigf
g(end-1, 1) = -.5*trace(dummy* (2/sigf*K *sigf)); % + 10*(sigf/signs)^9/maxSNR^10/signs *sigf;
% dK / dsign
g(end, 1) =  -.5*trace(2*signs*dummy*signs);  %- 10*(sigf/signs)^9/maxSNR^10*sigf/signs^2 *signs;
% dK / dwi
for i = 1:d
	Wloc = zeros(d); Wloc(i, i) = w(i)^-3;
	g(i, 1) = -.5*trace(dummy*((K.*maha(input, input, Wloc))* w(i) )); %+ 10*(1/std(input(:, i))/100)^10 * w(i)^9 *w(i);
end

if verb
	f, g
end
