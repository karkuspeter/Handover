function hypKernel = medianTrick(samples, quant)
if nargin < 2
    quant = 0.5;
end

for i = 1:size(samples, 2)
   
    dist = maha(samples(:, i), samples(:, i));
    hypKernel(i) = .5*median(quantile(dist.^.5, quant));
    
    
end