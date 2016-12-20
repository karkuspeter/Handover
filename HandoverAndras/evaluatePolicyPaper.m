
close all, clear all
addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

warning('off')

userData;

load('HandoverLearningOrientation_Peter.mat')

figure,
plot(data.meanR_initMean *9/4 + 5.5)
hold on, plot(data.meanR_learntMean *9/4 + 5.5)
legend('E[R] init', 'E[R] learnt' ) 

