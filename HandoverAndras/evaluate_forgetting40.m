clear all, close all,

addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

warning('off')

ridge = 1e-4;

fixedActivation = 0.2;

userData;
figure
for j =1:length(user.names)
    load(['HandoverLearningOrientation_', user.names{j}, '.mat'])
    
    kldiv(j) = data.KLdiv(end);
    diffReported(j) = mean(user.finalScores(j, :))-mean(user.initScores(j, :));
    dummy = corrcoef([data.meanR', data.meanR_init']);
    corr(j) = dummy(1, 2);
    dpol(j) = data.diffPolicy(end);
    dpolmean(j) = data.diffPolicyMean(end);
    meanKldiv(j) = mean(data.KLdiv);
    diffRModel(j) = data.meanR(end) - data.meanR_init(end);
    
    if diffReported(j) > 0
        marker{j} = 'o';
    elseif and(diffReported(j) > -.5, diffReported(j) < 0)
        marker{j} = '*';
    else
        marker{j} = 'x';
    end
    
    subplot(5, 2, 1)
     plot(diffReported(j),kldiv(j),  marker{j}), hold on
     ylabel('kl div'), xlabel('diff R reported')
     
     subplot(5, 2, 3)
     plot(diffReported(j), corr(j), marker{j}), hold on
     ylabel('corr policy mean'), xlabel('diff R reported')
     
     subplot(5, 2, 5)
     plot(diffReported(j), dpol(j), marker{j}), hold on
     ylabel('diff final policy R'), xlabel('diff R reported')
     
     subplot(5, 2, 7)
     plot(diffReported(j), dpolmean(j), marker{j}), hold on
     ylabel('diff final policy mean R'), xlabel('diff R reported')
     
     subplot(5, 2, 9)
     plot(diffReported(j),meanKldiv(j),  marker{j}), hold on
     ylabel('mean kl div'), xlabel('diff R reported')
     
      subplot(5, 2, 2)
     plot(diffRModel(j),kldiv(j),  marker{j}), hold on
     ylabel('kl div'), xlabel('diff R model')
     
     subplot(5, 2, 4)
     plot(diffRModel(j), corr(j), marker{j}), hold on
     ylabel('corr policy mean'), xlabel('diff R model')
     
     subplot(5, 2, 6)
     plot(diffRModel(j), dpol(j), marker{j}), hold on
     ylabel('diff final policy R'), xlabel('diff R model')
     
     subplot(5, 2, 8)
     plot(diffRModel(j), dpolmean(j), marker{j}), hold on
     ylabel('diff final policy mean R'), xlabel('diff R model')
     
     subplot(5, 2, 10)
     plot(diffRModel(j),meanKldiv(j),  marker{j}), hold on
     ylabel('mean kl div'), xlabel('diff R model')
    
end

