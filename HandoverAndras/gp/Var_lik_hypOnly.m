function [f, df] = Var_lik_hypOnly(w, px, x, y, ixPseudo, useGPU, useSinglePrecision)
%
% Calculates the bound and its derivates w.r.t. hyperparameters
%	(for details see 'Variational Learning of Inducing Variables in Sparse
%	 Gaussian Processes' by M.K. Titsias (AISTATS, 2009)
%	Note: - we assume pseudo inputs are given and is the same for all output dimensions
%		  - the algorithm calculates the gradients in 'one shot' for each output dimensions 
%			as they are independent given the pseudo inputs and w
%	      - we are using squared exponential kernel
%
% Inputs:
%	w: log hyperparameters (as in SPGP: b, c, sig) [(d+2)*e, 1]
%	px: pseudo inputs [m, d]
%	x: training inputs [n, d]
%	y: training targets [n, e]
%	ixPseudo: index that shows which samples are pseudo [n, 1] (logical)
%	useGPU: compute on GPU {[0], 1}
%	useSinglePrecision: uses single instead of double precision {[0], 1}
%
% Outputs:
%	f: neg lower bound ~ upper bound to minimize [1]
%	df: the gradients [(d+2)*e, 1]

usePaperVersion = 0;

if nargin < 5
	useGPU = 0;
	useSinglePrecision = 0;
end

[n, e] = size(y);
d = length(w)/e - 2;
m = size(px, 1);

hyp = reshape(w, d+2, e);
B = exp(hyp(1:d, :));
C = exp(hyp(d+1, :));
Sig = exp(hyp(d+2, :));

classType = 'double';
if useSinglePrecision
	classType = 'single';
end

Kmm = zeros(m, m, classType);
Knm = zeros(n, m, classType);
Qstar = zeros(n, n, classType);
Qnn = zeros(n, n, classType);
f = 0;
df = zeros(d+2, e, classType);

eyen = eye(n, classType);
eyem = eye(m, classType);
if useGPU
	Kmm = gpuArray(Kmm);
	Knm = gpuArray(Knm);
	Qstar = gpuArray(Qstar);
	Qnn = gpuArray(Qnn);
	px = gpuArray(px);
	x = gpuArray(x);
	y = gpuArray(y);
	B = gpuArray(B);
	C = gpuArray(C);
	Sig = gpuArray(Sig);
	eyen = parallel.gpu.GPUArray.eye(n, classType);
	eyem = parallel.gpu.GPUArray.eye(m, classType);
	f = gpuArray(f);
	df = gpuArray(df);
end


for i = 1:e
	Kmm = C(i)*exp(-.5*maha(px, px, diag(B(:, i)))) + 1e-6*eyem;
	Lmm = chol(Kmm);
	Knm = C(i)*exp(-.5*maha(x , px, diag(B(:, i))));
	Qnn1 = Knm/Lmm;
	Qnn = Qnn1*Qnn1';

	Qstar =  Qnn + eyen * Sig(i) + 1e-6*eyen;
	Lqs = chol(Qstar);
	iQstar = eyen/Lqs/Lqs';
	%iQstar = inv(single(Qstar));
	if usePaperVersion == 1
		QnnTruncated = Qnn(~ixPseudo, ~ixPseudo);
		diagKhat = ones(n-m, 1, classType)*C(i) - diag(QnnTruncated);
	else
		QnnTruncated = Qnn;
		diagKhat = ones(n, 1, classType)*C(i) - diag(QnnTruncated);
		ixPseudo = zeros(n, 1);
	end

	logDet = 2*sum(log(diag(Lqs)));

	f = f - n/2*log(2*pi) - .5*logDet - .5*y(:, i)'*iQstar*y(:, i) - ...
		.5/Sig(i)*sum(diagKhat);

	% dsig
	df(d+2, i) = -.5 * trace(iQstar) ...
				 +.5 * y(:, i)'*iQstar*iQstar*y(:, i) ...
				 +.5 / Sig(i) / Sig(i) * sum(diagKhat);
	% dlogsig
	df(d+2, i) = df(d+2, i) * Sig(i);

	% dc
	df(d+1, i) = -.5*sum(sum(bsxfun(@times, iQstar', Qnn)))/C(i) ...%-.5 * trace(iQstar*Qnn/C(i)) ...
				 +.5 * y(:, i)'*iQstar*Qnn*iQstar*y(:, i) / C(i) ...
				 -.5 / Sig(i) * sum(ones(size(QnnTruncated, 1), 1) - diag(QnnTruncated)/C(i));
	% dlogc
	df(d+1, i) = df(d+1, i) * C(i);
	
	% db
	dummy1 = Knm/Kmm;
	for j = 1:d
		dummy2 = -.5*Kmm.*(bsxfun(@minus, px(:, j), px(:, j)').^2);
		dummy3 = -.5*Knm.*(bsxfun(@minus, x (:, j), px(:, j)').^2);
		QnnDbi = dummy3*dummy1' - dummy1*dummy2*dummy1' + dummy1*dummy3';
		df(j, i) =  -.5 * sum(sum(bsxfun(@times, iQstar', QnnDbi))) ...
					+.5*y(:, i)'*iQstar*QnnDbi*iQstar*y(:, i) ...
					-.5/Sig(i) * sum(-diag(QnnDbi(~ixPseudo, ~ixPseudo)));
		df(j, i) = df(j, i) * B(j, i);
	end

end

df = -reshape(df, e*(d+2), 1);
f = -f;

if useGPU
	f = gather(f);
	df = gather(df);
end




