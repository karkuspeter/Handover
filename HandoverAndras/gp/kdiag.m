function Kd = kdiag(x,hyp);

c = exp(hyp(end-1));
% Kd = repmat(c,size(x,1),1);
Kd = c*ones(size(x, 1), 1);