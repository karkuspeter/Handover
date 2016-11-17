clear all, close all

data.fileName = 'HandoverLearning_test';
data.epsilon = 0.7;
data.initSamples = 29;
data.updateSamples = 5;
data.numHyper = 7;
data.samples = [];
data.policyMean = [.7 1.1 1.1 0 -150];
data.policyCov = {diag([.1,.3,.3,100,50].^2)};
data.policyStd = [.1,.3,.3,100,50];
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
    
