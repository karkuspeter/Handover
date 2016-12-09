clear all, close all

clear all, close all,

userData;

numAbsFeedback = [];

for i = 1:length(user.names)
    
    load(['HandoverLearningOrientation_', user.names{i}, '.mat']);
    
    numAbsFeedback(i) = length(data.absFeedback);
    
end

mean(numAbsFeedback)
std(numAbsFeedback)
