clear all, close all

data.fileName = 'HandoverLearning_test';
data.epsilon = 0.5;
data.initSamples = 20;
data.updateSamples = 5;
data.samples = [];
% data.policyMean = [.7 1.1 1.1 0 -150];
% data.policyMean =  [0.5,.5, .5, -300, -300];
data.policyMean =  [.8, 1.8, 1.8, -150, -150];
data.policyCov = {diag([.01,.01,.01,.01,.01].^2)};
data.prefFeedback = [  ];
data.absFeedback = [  ];
data.failedExperiments = [];
data.hyp = [];
data.polparMinLimit = [0.5,.5, .5, -300, -300];
data.polparMaxLimit = [.95, 1.8, 1.8, 300, 0];

if exist(['./', data.fileName, '.mat'], 'file')
    error([data.fileName, '.mat already exists!'])
else
    save(['./',data.fileName, '.mat'], 'data');
end
    
