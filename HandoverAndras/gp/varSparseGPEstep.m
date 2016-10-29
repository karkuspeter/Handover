function Fv = varSparseGPEstep(w, px, x, y, ixPseudo)
%
% Calculates the lower bound when adding a non-pseudo input to the pseudo database
%
% Inputs:
%	hyp: log hyperparameters as in SPGP (b, c, sig) [d+2, e];
%	px: pseudo inputs [m, d]
%	x: full training inputs [n, d]
%	y: full training targets [n, e]
%	ixPseudo: the index of pseudo Inputs [n, 1]
%
% Output
%	Fv: the new lower bound after including each non-pseudo trining input to the pseudo database [n-m, 1]

usePaperVersion = 0;

[n, d] = size(x);
m = size(px, 1);
e = size(y, 2);

hyp = reshape(w, d+2, e);
B = exp(hyp(1:d, :));
C = exp(hyp(d+1, :));
Sig = exp(hyp(d+2, :));

classType = 'double';
%if useSinglePrecision
%	classType = 'single';
%end

Kmm = zeros(m, m, e, classType);
Knm = zeros(n, m, e, classType);
Qstar = zeros(n, n, classType);
Qnn = zeros(n, n, classType);
diagKhat = zeros(n-m-1, 1, classType);
f = 0;
df = zeros(d+2, e, classType);

Fv = zeros(n-m, 1, classType);
ixNonPseudo = find(ixPseudo == 0);

for i = 1:e
	Knn = C(i)*exp(-.5*maha(x , x, diag(B(:, i)))) + 1e-6*eye(n, classType);
	Kmm = Knn(ixPseudo, ixPseudo);
	Lmm = chol(Kmm);
	iKmm = eye(m, classType)/Lmm/Lmm';
	Knm = Knn(:, ixPseudo);
	%Kmm = C(i)*exp(-.5*maha(px, px, diag(B(:, i)))) + 1e-6*eye(m, classType);
	%Lmm = chol(Kmm);
	%Knm = C(i)*exp(-.5*maha(x , px, diag(B(:, i))));
	for j = 1:(n-m)
		ixNew = ixNonPseudo(j);
		colNew = Knn(ixPseudo, ixNew);
		diagNew = C(i);
		iKmmCurr = covInverseUpdate(iKmm, colNew, diagNew);
		KnmCurr = [Knm, Knn(:, ixNew)];
		Qnn = KnmCurr*iKmmCurr*KnmCurr';
		ixPseudoCurr = ixPseudo;
		ixPseudoCurr(ixNew) = 1 == 1;

		Qstar =  Qnn + eye(n, classType) * Sig(i);
		Lqs = chol(Qstar);
		if usePaperVersion == 1
			QnnTruncated = Qnn(~ixPseudoCurr, ~ixPseudoCurr);
			diagKhat = ones(n-m-1, 1, classType)*C(i) - diag(QnnTruncated);
		else
			QnnTruncated = Qnn;
			diagKhat = ones(n, 1, classType)*C(i) - diag(QnnTruncated);
		end

		logDet = 2*sum(log(diag(Lqs)));

		Fv(j) = Fv(j) - n/2*log(2*pi) - .5*logDet - .5*y(:, i)'/Lqs/Lqs'*y(:, i) - ...
		.5/Sig(i)*sum(diagKhat);
	end
end


