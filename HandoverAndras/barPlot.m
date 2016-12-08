clear all, close all,

userData;

barData = [];
errorBarData = [];

for i = 1:length(user.names)
    
    load(['HandoverLearningOrientation_', user.names{i}, '.mat']);
    learntMean = data.meanR_learntMean(end)* 9/4 + 5.5;
    initMean = data.meanR_initMean(end)*9/4 + 5.5;
    barData(i, :) = [learntMean, mean(user.finalScores(i, :)), initMean, mean(user.initScores(i, :))];
    errorBarData(i, 1) = data.stdR_learntMean(end)*2*(9/4);
    errorBarData(i, 2) = std(user.finalScores(i, :));
    errorBarData(i, 3) = data.stdR_initMean(end)*2*(9/4);
    errorBarData(i, 4) = std(user.initScores(i, :));
    
    
end

figure, bar(barData)
legend('LearntRewardMean', 'LearntReportedMean', 'InitRewardMean', 'ReportedInitMean')
xlabel('#subject')
ylabel('Reward')
title('Reward model prediction vs reported rewards')

figure, barwitherr(errorBarData, barData)
legend('LearntRewardMean', 'LearntReportedMean', 'InitRewardMean', 'ReportedInitMean')
xlabel('#subject')
ylabel('Reward')
title('Reward model prediction vs reported rewards')