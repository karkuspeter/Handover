function [f, dfX] = Var_lik(W, X, Xu, y, useGPU, useSinglePrecision)
%[f, dfX] = Var_lik_Hyp(W, X, y, m)
%
%Description:  Computes the negative lower bound value and
%              the gradient over pseudo-inputs and hyperparameters
%              for the squared exponential kernel with varied
%              lengthscales
%
% Inputs:
%        m: number of pseudo points
%        X: inputs
%        Xu: pseudo inputs
%        y: outputs
%        W: stores the hyperparameters
%           logtheta stores the hyperparameters exactly in the same
%           format as the software of R&W (2006)
%
% Outputs:
%         f: minus the lower bound value
%         dfX: vectors with the derivatives stored as in W vector

% number of examples and dimension of input space
m = size(Xu, 1);
if useSinglePrecision
	classType = 'single';
else
	classType = 'double';
end

if useGPU
    W = gpuArray(W);
    X = gpuArray(X);
    y = gpuArray(y);
    classType = 'double';
    if useSinglePrecision
        W = single(W);
        X = single(X);
        y = single(y);
        classType = 'single';
    end
	eyem = parallel.gpu.GPUArray.eye(m, classType);
else
	eyem = eye(m, classType);
end

[n, D] = size(X);
jitter = 1e-6;

% extract pseudo inptus and hyperparameters
logtheta = W(1:D+2);
sigma2n = exp(2*logtheta(D+2));
sigmaf = exp(2*logtheta(D+1));

X = X ./ repmat(exp(logtheta(1:D))',n,1);
Xu = Xu ./ repmat(exp(logtheta(1:D))',m,1);

% covariances
Kmm = Xu*Xu';
Kmm = repmat(diag(Kmm),1,m) + repmat(diag(Kmm)',m,1) - 2*Kmm;
Kmm = sigmaf*exp(-0.5*Kmm);
Kmm = Kmm + (jitter*sigmaf)*eyem;
Knm = -2*Xu*X' + repmat(sum(X.*X,2)',m,1) + repmat(sum(Xu.*Xu,2),1,n);
Knm = sigmaf*exp(-0.5*Knm');

try
  Lkmm = chol(Kmm);
catch
  f = Inf;
  dfX = zeros(size(W));
  return
end
Cnm1 = Knm/Lkmm;
Cmnmn = Cnm1'*Cnm1; %Qnn
Lm = chol(sigma2n*eyem + Cmnmn);
invLm = Lm\eyem;
Pnm1 = Cnm1/Lm;

Pmnmn = (Lm'\Cmnmn)';

bet = Pnm1'*y;

LmLkmm = Lm*Lkmm;
wm1 = (LmLkmm)\bet;

invQt1 = y - Pnm1*bet;

logdetQ = (n-m)*2*logtheta(D+2) + 2*sum(log(diag(Lm)));
f = 0.5*logdetQ +  (0.5/sigma2n)*(y'*y - bet'*bet)  + 0.5*n*log(2*pi);
% add the trace term
TrK = + (0.5/sigma2n)*(n*sigmaf - sum(diag(Cmnmn)));
f = f + TrK;

if nargout>1
  % compute derivatives
    if useGPU
		df = parallel.gpu.GPUArray.zeros(D+2,1, classType);
	else
		df = zeros(D+2,1);
	end

  % precomputations
  aux = sum(sum(invLm.*Pmnmn));
  Pmnmn = (Lkmm\Cmnmn)';
  BB1 = Lkmm\Pmnmn;
  BB1 = - BB1/sigma2n - wm1*wm1';
  BB1 = BB1 + LmLkmm\(Lm'\Pmnmn);
  BB1 = Kmm.*BB1;

  Cnm1 = (Lkmm\Cnm1')';
  Cnm1 = Cnm1 - sigma2n*(LmLkmm\Pnm1')';
  Pnm1 = repmat(invQt1,1,m).*repmat(wm1',n,1);

  Cnm1 = (Cnm1 + Pnm1).*Knm;
  Pnm1 = Pnm1.*Knm;
  %
  for d = 1:D
    Knm = -maha(X(:,d),Xu(:,d));
    Kmm = -maha(Xu(:,d),Xu(:,d));

    df(d) = sum(sum(Knm.*Cnm1))/sigma2n + 0.5*sum(sum(Kmm.*BB1));
  end
  %
  df(D+1) = -sum(sum(Pnm1))/sigma2n+aux;
  df(D+2) = (n-aux) - (invQt1'*invQt1)/sigma2n;

  % trace term contribution to the sigmaf and sigma2n
  df(D+1) = df(D+1) + 2*TrK;
  df(D+2) = df(D+2) - 2*TrK;

  dfX = df;
end

if useGPU
	dfX = gather(dfX);
	f = gather(f);
end

