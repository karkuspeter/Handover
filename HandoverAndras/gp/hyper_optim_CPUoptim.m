function [f, g] = hyper_optim_CPUoptim(loghyp, input, target, verb)


% loghyp: the log-hyperparameters log(w, sigf, sign)
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
L = chol(Ky);

alpha =  L\(L'\target);

ldKy = 2*sum(log(diag(L)));
f = -(-.5*target'*alpha-.5*ldKy - n/2*log(2*pi)) ;

dummy = alpha*alpha';
dummy = dummy -eye(n)/L/L';

g = zeros(length(loghyp), 1);

k = 2; th = 100;
SNR = sigf/signs;
SNRpenalty = -log(factorial(k)) -k*log(th) + (k-1)*log(SNR) - SNR/th;
scale = -length(target)/5;

if SNR < 500
    scale = 0;
end

f = f + scale*SNRpenalty;

dSNRpenalty = scale*(-1/th + (k-1)/SNR)*SNR/2;

dsign_snr = -dSNRpenalty;
dsigf_snr = dSNRpenalty;


% dK / dsigf
p1 = dummy;
p1 = p1*(2/sigf*K *sigf);

g(end-1, 1) = -.5*trace(p1) + dsigf_snr; 
% dK / dsign
g(end, 1) =  -.5*trace(2*signs*dummy*signs)  + dsign_snr; 
% dK / dwi

for i = 1:d
	Wloc = zeros(d); Wloc(i, i) = w(i)^-3;
    
    p1 = maha(input, input, Wloc);
    p1 = (K.* p1)* w(i);
    p1 = dummy*p1;
	g(i, 1) = -.5*trace(p1); 
end

if verb
	f, g
end

end

function [f, df] = gamma_func(x, k, th)
%
% The gamma function and its derivate w.r.t. x
%
% input:
%   k, th: gamma function parameters
%   x: x axis parameters (x > 0) [n, 1]
% output:
%   f: the gamma functin [n, 1]
%   df: the derivate of gamma func [n, 1]

f = 1/factorial(k)/th^k * x.^(k-1) .* exp(-x./th);
df = 1/factorial(k)/th^k .* exp(-x./th) .* (...
    (k-1)*x.^(k-2) - 1/th* x.^(k-1));
end
