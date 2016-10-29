clear all, close all
addpath ~/svnprojects/ClassSystem/Helpers/shadedErrorBar/
load Shandover

for i = 1:20
    jani = length(S.rewMean{i});
    if jani < 6
        S.rewMean{i} = [S.rewMean{i}, S.rewMean{i}(end) * ones(6-jani, 1)'];
        S.rewStd{i} = [S.rewStd{i}, S.rewStd{i}(end) * ones(6-jani, 1)'];
        S.rewMeanAbs{i} = [S.rewMeanAbs{i}, S.rewMeanAbs{i}(end) * ones(6-jani, 1)'];
        for j = jani+1:6
            S.aA{i}(:, :,j)  = S.aA{i}(:, :, jani);
        end
    end
end
        
R = cell2mat(S.rewMean');
% Rabs = cell2mat(S.rewMeanAbs');


fs = 24;
ms = 16;

figure, 
shadedErrorBar([40 50 60 70 80 90], mean(R, 1), std(R))
xlabel('Robot Experiments', 'FontSize', fs)
ylabel('Reward', 'FontSize', fs)
title('Human Robot Handover Learning', 'FontSize', fs)

set(gca,'FontSize',fs)

% 
% 
% 
% 
% jani = zeros(7, 2, 20);
% for i = 1:20
%    jani(:, :, i) = S.aA{i}(:, :, end);
% end