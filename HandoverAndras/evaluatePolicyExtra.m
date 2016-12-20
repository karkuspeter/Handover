close all, clear all
addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

warning('off')

ridge = 1e-4;

fixedActivation = 0.2;

userData;

for j =1:length(user.names)
    load(['HandoverLearningOrientation_', user.names{j}, '.mat'])
    
    clf,
    plot(data.meanR), hold on
    plot(data.meanR_init),
    plot(data.meanR_initMean),
    plot(data.meanR_learntMean)
    pause
end