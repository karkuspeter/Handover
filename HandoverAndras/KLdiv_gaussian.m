function kldiv = KLdiv_gaussian(mu1, cov1, mu2, cov2)

mu1 = mu1(:);
mu2 = mu2(:);   
k = length(mu1);
kldiv = .5 * (trace(cov2\cov1) + (mu2-mu1)'/cov2*(mu2-mu1) -k + log(det(cov2)/det(cov1)));