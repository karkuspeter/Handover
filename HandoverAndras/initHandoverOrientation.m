clear all, close all

data.fileName = 'HandoverLearningOrientation_test';
data.epsilon = 0.7;
data.initSamples = 29;
data.updateSamples = 5;
data.numHyper = 6;
data.samples = [];
data.policyMean = [.7 1.1 0, 0];
% data.policyMean = [.7 1.1 0, -20/180*pi];
% data.policyCov = {diag([.1,.3,50, .05].^2)};
data.policyCov = {diag([.1,.3,50, 30/180*pi].^2)};
data.policyStd = [.1,.3,50, 30/180*pi];
data.prefFeedback = [  ];
data.absFeedback = [  ];
data.failedExperiments = [];
data.hyp = [];
data.polparMinLimit = [0.5, .5, -300,  -90/180*pi ];
data.polparMaxLimit = [.95, 1.8, 300, 90/180*pi];
data.meanR = [];
data.stdR = [];

if exist(['./', data.fileName, '.mat'], 'file')
    error([data.fileName, '.mat already exists!'])
else
    save(['./',data.fileName, '.mat'], 'data');
end
    
