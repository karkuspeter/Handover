clear all, close all
addpath('Batch 11')
addpath('Batch 12')
addpath('Batch 13')
addpath('Batch 14')
num = 196;
[fPrior, prefs, context, samples] = humanFeedbackRobotHandover();
disp(['context: ', num2str(context(num-16))])
jani = dlmread([num2str(num), '.txt']);

posX = jani(:, 15);
posY = jani(:, 16);
posZ = jani(:, 17);

ixStart = find(and(and(posX > 750, posY > -250), posZ < 100));
ixStart = 1;
ixEnd = size(jani, 1);
ixOK = ixStart:ixEnd;

time = jani(ixOK, 1) + jani(ixOK, 2)/1e9;
time = time - time(1);
forceMagn = smooth(sqrt(jani(ixOK, 5).^2 + .25 * jani(ixOK, 7).^2));
torqueMagn = sqrt(jani(ixOK, 9).^2+jani(ixOK, 11).^2+jani(ixOK, 13).^2);
jerkMagn = diff(forceMagn)/.02;
posX = jani(ixOK, 15);
posY = jani(ixOK, 16);
posZ = jani(ixOK, 17);
posMagn = sqrt(posX.^2 + posY.^2 + posZ.^2);
velMagn = diff(smooth(posMagn))/.02;
fs =18;

figure,
subplot(2,2,1)
plot(time, forceMagn);ylabel('Abs. Force [N]', 'FontSize', fs)

set(gca,'FontSize',fs)
subplot(2,2,2)
plot(time(1:end-1), jerkMagn);ylabel('Abs. Jerk [N/s]', 'FontSize', fs)
set(gca,'FontSize',fs)
subplot(2,2,3)
plot(time, smooth(jani(ixOK, 3))), ylabel('Force_X [N]', 'FontSize', fs)
set(gca,'FontSize',fs)
subplot(2,2,4)
plot(time(1:end-1), diff(smooth(jani(ixOK, 3)))/.02), ylabel('Jerk_X [N/s]', 'FontSize', fs)
set(gca,'FontSize',fs)

