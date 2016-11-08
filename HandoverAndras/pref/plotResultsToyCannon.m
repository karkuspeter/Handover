clear all, close all
figure,
addpath('~/Downloads')

% random 3 feedback + 20% absolute feedback

fs = 18;
ms = 12;

colors = {'b-', 'r-', 'k-', 'b--', 'r--', 'k--'};
symbols = ['o', '<', '+', 's', '>', 'x'];

load PrefToyCannonSimple3_Consec0_AbsFreq5_InitNewSamples40_5_eps75_siga0
numUpdates = size(R.mean, 2);
x = 40 + [0:size(R.mean, 2)-1] * 5;
N = 3;
xN = x(1:N:end);

figure(1)
plot(xN, smooth(mean(R.mean(:, 1:N:end)))-4, [colors{1}, symbols(1)], 'MarkerSize', ms, 'MarkerEdgeColor', colors{1}(1), 'MarkerFaceColor', 'none'), hold on,
figure(2)
plot(xN, smooth(std(R.mean(:, 1:N:end))), [colors{1}, symbols(1)], 'MarkerSize', ms, 'MarkerEdgeColor', colors{1}(1), 'MarkerFaceColor', 'none'), hold on

% load RepsToyCannon_NewSamp40_eps75.mat
% figure(1)
% plot(xN, smooth(mean(R.mean, 1))-4, [colors{2}, symbols(2)], 'MarkerSize', ms, 'MarkerEdgeColor', colors{2}(1), 'MarkerFaceColor', 'none'), hold on,
% figure(2)
% plot(xN, smooth(std(R.mean)), [colors{2}, symbols(2)], 'MarkerSize', ms, 'MarkerEdgeColor', colors{2}(1), 'MarkerFaceColor', 'none'), hold on
% keyboard
load PrefToyCannonSimple1_Consec1_AbsFreq5_InitNewSamples40_5_eps75_siga0

figure(1)
plot(xN, smooth(mean(R.mean(:, 1:N:end)))-4, [colors{2}, symbols(2)], 'MarkerSize', ms, 'MarkerEdgeColor', colors{2}(1), 'MarkerFaceColor', 'none'), hold on,
figure(2)
plot(xN, smooth(std(R.mean(:, 1:N:end))), [colors{2}, symbols(2)], 'MarkerSize', ms, 'MarkerEdgeColor', colors{2}(1), 'MarkerFaceColor', 'none'), hold on

% consecutive feedback + 20% absolute feedback
load PrefToyCannonSimple1_Consec0_AbsFreq5_InitNewSamples40_5_eps75_siga0

figure(1)
plot(xN, smooth(mean(R.mean(:, 1:N:end)))-4, [colors{3}, symbols(3)], 'MarkerSize', ms, 'MarkerEdgeColor', colors{3}(1), 'MarkerFaceColor', 'none'), hold on,
figure(2)
plot(xN, smooth(std(R.mean(:, 1:N:end))), [colors{3}, symbols(3)], 'MarkerSize', ms, 'MarkerEdgeColor', colors{3}(1), 'MarkerFaceColor', 'none'), hold on

load PrefToyCannonSimple3_Consec0_AbsFreq10000_InitNewSamples40_5_eps75_siga0

figure(1)
plot(xN, smooth(mean(R.mean(:, 1:N:end)))-4, [colors{4}, symbols(4)], 'MarkerSize', ms, 'MarkerEdgeColor', colors{4}(1), 'MarkerFaceColor', 'none'), hold on,
figure(2)
plot(xN, smooth(std(R.mean(:, 1:N:end))), [colors{4}, symbols(4)], 'MarkerSize', ms, 'MarkerEdgeColor', colors{4}(1), 'MarkerFaceColor', 'none'), hold on

load PrefToyCannonSimple1_Consec1_AbsFreq10000_InitNewSamples40_5_eps75_siga0

figure(1)
plot(xN, smooth(mean(R.mean(:, 1:N:end)))-4, [colors{5}, symbols(5)], 'MarkerSize', ms, 'MarkerEdgeColor', colors{5}(1), 'MarkerFaceColor', 'none'), hold on,
figure(2)
plot(xN, smooth(std(R.mean(:, 1:N:end))), [colors{5}, symbols(5)], 'MarkerSize', ms, 'MarkerEdgeColor', colors{5}(1), 'MarkerFaceColor', 'none'), hold on

% consecutive feedback + 20% absolute feedback
load PrefToyCannonSimple1_Consec0_AbsFreq10000_InitNewSamples40_5_eps75_siga0

figure(1)
plot(xN, smooth(mean(R.mean(:, 1:N:end)))-4, [colors{6}, symbols(6)], 'MarkerSize', ms, 'MarkerEdgeColor', colors{6}(1), 'MarkerFaceColor', 'none'), hold on,
figure(2)
plot(xN, smooth(std(R.mean(:, 1:N:end))), [colors{6}, symbols(6)], 'MarkerSize', ms, 'MarkerEdgeColor', colors{6}(1), 'MarkerFaceColor', 'none'), hold on

figure(1), 
h = legend('Random 3 pref, 1 in 5 absolute', 'Consequent pref, 1 in 5 absolute', 'Random 1 pref, 1 in 5 absolute', 'Random 3 pref', 'Consequent pref', 'Random 1 pref');
set(h, 'Box', 'off')
set(h, 'FontSize', fs)
set(gca,'FontSize',fs)
title('Toy Cannon task', 'FontSize', fs)

xlabel('Sample evaluations','FontSize', fs)
ylabel('Expected reward','FontSize', fs)
axis([38 190, -6.0 0])

figure(2), 
% h = legend('Random 3 pref, 1/5 absolute', 'Consecuent pref, 1/5 absolute', 'Random 1 pref, 1/5 absolute', 'Random 3 pref', 'Consecuent pref', 'Random 1 pref');
% set(h, 'FontSize', fs)
set(gca,'FontSize',fs)
axis([38 190, 0, 2])
title('Toy Cannon task', 'FontSize', fs)
xlabel('Sample evaluations','FontSize', fs)
ylabel('STD of Expected reward','FontSize', fs)

load PrefToyCannonSimple1_Consec1_AbsFreq5_InitNewSamples40_5_eps75_siga0

figure(3)
plot(xN, smooth(mean(R.mean(:, 1:N:end)))-4, ['k-<'], 'MarkerSize', ms, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none'), hold on,
figure(4)
plot(xN, smooth(std(R.mean(:, 1:N:end))), ['k-<'], 'MarkerSize', ms, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'none'), hold on,
load  PrefToyCannonSimpleGP_AbsFreq1_InitNewSamples40_5_eps75_siga0.mat             
figure(3)
plot(xN, smooth(mean(R.mean(:, 1:N:end)))-4,  ['m-d'], 'MarkerSize', ms, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'none'), hold on,
figure(4)
plot(xN, smooth(std(R.mean(:, 1:N:end))), ['m-d'], 'MarkerSize', ms, 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'none'), hold on,
size(R.mean, 1)

load  PrefToyCannonSimpleGP_AbsFreq3_InitNewSamples40_5_eps75_siga0.mat             
figure(3)
plot(xN, smooth(mean(R.mean(:, 1:N:end)))-4, ['b-o'], 'MarkerSize', ms, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'none'), hold on,
figure(4)
plot(xN, smooth(std(R.mean(:, 1:N:end))),  ['b-o'], 'MarkerSize', ms, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'none'), hold on,
size(R.mean, 1)
load PrefToyCannonSimpleGP_AbsFreq5_InitNewSamples40_5_eps75_siga0.mat
figure(3)
plot(xN, smooth(mean(R.mean(:, 1:N:end)))-4, ['r-s'], 'MarkerSize', ms, 'MarkerEdgeColor','r', 'MarkerFaceColor', 'none'), hold on,
figure(4)
plot(xN, smooth(std(R.mean(:, 1:N:end))),  ['r-s'], 'MarkerSize', ms, 'MarkerEdgeColor','r', 'MarkerFaceColor', 'none'), hold on,


figure(3)
h = legend('Consequent pref, 1 in 5 absolute', '1 in 1 absolute', '1 in 3 absolute', '1 in 5 absolute');
set(h, 'Box', 'off')
set(h, 'FontSize', fs)
set(gca,'FontSize',fs)
ylabel('Expected reward', 'FontSize', fs)
xlabel('Sample evaluations', 'FontSize', fs)
title('Toy Cannon task', 'FontSize', fs)
axis([38 190, -4.5 0])
figure(4)
% h = legend('Consecuent + 1/5 absolute', 'GP-1', 'GP-3', 'GP-5');
% set(h, 'FontSize', fs)
set(gca,'FontSize',fs)
ylabel('STD of Expected reward', 'FontSize', fs)
xlabel('Sample evaluations', 'FontSize', fs)
title('Toy Cannon task', 'FontSize', fs)
axis([38 190, 0, 2])


% 
% [colors(i),'-', symbols(i)], 'MarkerSize', ms, 'MarkerEdgeColor', colors(i), 'MarkerFaceColor', 'none');
% 
% h = legend('GPRL', 'GPTD', 'L22-KLSTD', 'VFP');
%     ylabel('Expected MSE', 'FontSize', fs)
%     title('State Chain', 'FontSize', fs)
% set(h, 'FontSize', fs)
% set(gca,'FontSize',fs)
% 
% ylabel('Std MSE', 'FontSize', fs)
% xlabel('Observed trajectories', 'FontSize', fs)
% set(gca,'FontSize',fs)