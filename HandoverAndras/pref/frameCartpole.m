clear all, close all

K =  [ -2
    -5
    -100
    -20];

dt = .02;
tend = 2;
forceStd = 1;
numTrajs = 5;

[r, xu] = simulateCartPole([.9 0 -1 0]', numTrajs, dt, tend, forceStd, K);
sum(sum(r))

figure, plot(r)
title('shitty controller')
figure, 
for i = 1:numTrajs
hold on, plot([xu(:, 1, i), xu(:, 3, i)])
end
title('shitty controller')

[r, xu] = simulateCartPole([.9 0 -1 0]', numTrajs, dt, tend, forceStd);
sum(sum(r))

figure, plot(r)
title('>>optimal<< controller')
figure, 
for i = 1:numTrajs
hold on, plot([xu(:, 1, i), xu(:, 3, i)])
end
title('>>optimal<< controller')

