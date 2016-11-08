clear all, close all

trajs = load('trajs');
trajs = trajs.trajs;

figure,
for i = 1:(size(trajs, 2)/12/2)
    ix_q = (i-1)*12+1:2:i*12-4;
    ix_dq = (i-1)*12+2:2:i*12-4;
    ix_u = (i-1)*12+9:i*12;
    subplot(3,1,1)
    plot(trajs(:, ix_q)), hold on
    subplot(3,1,2)
    plot(trajs(:, ix_dq)), hold on,
    subplot(3,1,3)
    plot(trajs(:, ix_u)), hold on
end