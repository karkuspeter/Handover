function hypKernel = medianTrick(samples, medianScale)
if nargin < 2
    medianScale = 0.1;
end

for i = 1:size(samples, 2)
   
    dist = maha(samples(:, i), samples(:, i));
    hypKernel(i) = medianScale*median(reshape(dist, [], 1));
    
    
end