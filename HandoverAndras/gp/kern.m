function K = kern(x1,x2,hyp)

[n1,dim] = size(x1); n2 = size(x2,1);
b = exp(hyp(1:end-2)); c = exp(hyp(end-1));

x1 = x1.*repmat(sqrt(b)',n1,1);
x2 = x2.*repmat(sqrt(b)',n2,1);
K = -2*x1*x2' + bsxfun(@plus, sum(x2.*x2,2)', sum(x1.*x1,2));
K = c*exp(-0.5*K);
% K = c*gauss_Knm(x1,x2,sqrt(2));