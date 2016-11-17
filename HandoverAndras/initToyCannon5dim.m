clear all, close all

data.fileName = 'ToyCannon5dim';
data.epsilon = 0.7;
data.initSamples = 30;
data.updateSamples = 5;
data.gradRestarts = 10;
data.numHyper = 7;
data.samples = [];
data.policyMean = [1, .1, .01, 1, .02];
data.policyCov = {diag(ones(1, 5)*.2).^2};
data.policyStd = ones(1, 5)*.2;
data.prefFeedback = [  ];
data.absFeedback = [];
data.failedExperiments = [];
data.hyp = [];
data.meanR = [];
data.stdR = [];

if exist(['./', data.fileName, '.mat'], 'file')
    error([data.fileName, '.mat already exists!'])
else
    save(['./',data.fileName, '.mat'], 'data');
end
    
