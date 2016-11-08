function logHypOpt = getFullGPModelGPUpso(input, target, n, steps)

[N, d] = size(input);
[N, e] = size(target);

wmin = std(input)/100;
wmax = std(input)*100;

for i = 1:e
    f = @(loghyp) hyper_optim_GPUoptim(loghyp, input, target(:, i), 0);
    
    range = log([wmin(:), wmax(:); std(target(:, i))/10, std(target(:, i))*50; std(target(:, i))/5000 std(target(:, i))/10]);
    
    logHypOpt(:, i) = pso(range, n, f, steps, [.5, .25, .25]);
end