function plotSampleTraj(nosamp, res, observed_sample_round)

osr = observed_sample_round;
len = res.len;
ix = ((nosamp-1)*len+1):(nosamp*len);

subplot(6,1,1)
plot(res.x(ix, 1:2:8), '--')
hold on, plot(res.traj_mean(ix, 1:2:8, osr));
ylabel('angle mean')

subplot(6,1,2)
plot(res.traj_std(ix, 1:2:8, osr));
ylabel('angle std')

subplot(6,1,3)
plot(res.x(ix, 2:2:8), '--')
hold on, plot(res.traj_mean(ix, 2:2:8, osr));
ylabel('angular vel mean')

subplot(6,1,4)
plot(res.traj_std(ix, 2:2:8, osr));
ylabel('angular vel std')

subplot(6,1,5)
plot(res.x(ix, 9:12), '--')
hold on, plot(res.traj_mean(ix, 9:12, osr));
ylabel('torque')

subplot(6,1,6)
plot(res.traj_std(ix, 9:12, osr));
ylabel('torque std')