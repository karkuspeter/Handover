function [M, S, V] = gp1(gpmodel, m, s)

% Compute joint predictions for the FITC sparse approximation to multiple GPs
% with uncertain inputs. If dynmodel.noise exists, individual noise contribu-
% tions are added. 
%
% dynmodel  dynamics model struct
%   hyp     log-hyper-parameters                                    [D+2 x  E ]
%   inputs  training inputs                                         [ n  x  D ]
%   target  training targets                                        [ n  x  E ]
%   induce  inducint inputs                                         [ np x  D ]
%   noise   (optional) individual noise variance terms              [ n  x  E ]
% m         mean of the test distribution                           [ D       ]
% s         covariance matrix of the test distribution              [ D  x  D ]
%
% M         mean of pred. distribution                              [ E       ]
% S         covariance of the pred. distribution                    [ E  x  E ]
% V         inv(s) times covariance between input and output        [ D  x  E ]
%
% Copyright (C) 2008-2012 by Marc Deisenroth & Carl Edward Rasmussen 2012-06-29

if ~isfield(gpmodel,'induce') || numel(gpmodel.induce)==0, 
    [M, S, V] = gp0(gpmodel, m, s); return; end

persistent iK iK2 beta oldX;
ridge = 1e-6;                        % jitter to make matrix better conditioned
[n, D] = size(gpmodel.inputs);    % number of examples and dimension of inputs
E = size(gpmodel.target,2);         % number of examples and number of outputs
X = gpmodel.hyp; input = gpmodel.inputs; target = gpmodel.target;

[np pD pE] = size(gpmodel.induce);     % number of pseudo inputs per dimension
pinput = gpmodel.induce;                                   % all pseudo inputs

if numel(X) ~= numel(oldX) || isempty(iK) || isempty(iK2) || ... % if necessary
              sum(any(X ~= oldX)) || numel(iK2) ~=E*np^2 || numel(iK) ~= n*np*E
  oldX = X;                                        % compute K, inv(K), inv(K2)
  iK = zeros(np,n,E); iK2 = zeros(np,np,E); beta = zeros(np,E);
    
  for i=1:E
    pinp = bsxfun(@rdivide,pinput(:,:,min(i,pE)),exp(X(1:D,i)'));
    inp = bsxfun(@rdivide,input,exp(X(1:D,i)'));
    Kmm = exp(2*X(D+1,i)-maha(pinp,pinp)/2) + ridge*eye(np);  % add small ridge
    Kmn = exp(2*X(D+1,i)-maha(pinp,inp)/2);
    L = chol(Kmm)';
    V = L\Kmn;                                             % inv(sqrt(Kmm))*Kmn
    if isfield(gpmodel,'noise')
      G = exp(2*X(D+1,i))-sum(V.^2)+gpmodel.noise(:,i)';
    else
      G = exp(2*X(D+1,i))-sum(V.^2);
    end
    G = sqrt(1+G/exp(2*X(D+2,i)));
    V = bsxfun(@rdivide,V,G);
    Am = chol(exp(2*X(D+2,i))*eye(np) + V*V')';
    At = L*Am;                                    % chol(sig*B) [thesis, p. 40]
    iAt = At\eye(np);
% The following is not an inverse matrix, but we'll treat it as such: multiply
% the targets from right and the cross-covariances left to get predictive mean.
    iK(:,:,i) = ((Am\(bsxfun(@rdivide,V,G)))'*iAt)';
    beta(:,i) = iK(:,:,i)*target(:,i);      
    iB = iAt'*iAt.*exp(2*X(D+2,i));              % inv(B), [Ed's thesis, p. 40]
    iK2(:,:,i) = Kmm\eye(np) - iB; % covariance matrix for predictive variances       
  end
  
%     feri = load('samplvsmatching');
%     iK2 = feri.s.vf;
%     beta = feri.s.mf';
end




k = zeros(np,E); M = zeros(E,1); V = zeros(D,E); S = zeros(E);       % allocate
inp = zeros(np,D,E);

for i=1:E     % compute predicted mean and inv(s) times input-output covariance
  inp(:,:,i) = bsxfun(@minus,pinput(:,:,min(i,pE)),m');
 
  L = diag(exp(-X(1:D,i)));
  in = inp(:,:,i)*L;
  B = L*s*L+eye(D); 
  
  t = in/B;
  l = exp(-sum(in.*t,2)/2); lb = l.*beta(:,i);
  tL = t*L;
  c = exp(2*X(D+1,i))/sqrt(det(B));
  
  M(i) = sum(lb)*c;                                            % predicted mean
  V(:,i) = tL'*lb*c;                     % inv(s) times input-output covariance
  k(:,i) = 2*X(D+1,i)-sum(in.*in,2)/2;
end

for i=1:E           % compute predictive covariance matrix, non-central moments
  ii = bsxfun(@rdivide,inp(:,:,i),exp(2*X(1:D,i)'));
  
  
  for j=1:i
    R = s*diag(exp(-2*X(1:D,i))+exp(-2*X(1:D,j)))+eye(D); t = 1./sqrt(det(R));
    ij = bsxfun(@rdivide,inp(:,:,j),exp(2*X(1:D,j)'));
    tmp = R\s/2;
    tmp2 = maha(ii,-ij,tmp);
%     L = exp(bsxfun(@plus,k(:,i),k(:,j)')+maha(ii,-ij,R\s/2));
    L = exp(bsxfun(@plus,k(:,i),k(:,j)')+tmp2);
    if i==j
      S(i,i) = t*(beta(:,i)'*L*beta(:,i) - sum(sum(iK2(:,:,i).*L)));
    else
      S(i,j) = beta(:,i)'*L*beta(:,j)*t; S(j,i) = S(i,j);
    end  
  end

  S(i,i) = S(i,i) + exp(2*X(D+1,i));
end

S = S - M*M';                                               % centralize moments
