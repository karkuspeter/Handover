clear all, close all

data.fileName = 'ToyCannon';
data.epsilon = 0.5;
data.initSamples = 20;
data.updateSamples = 5;
data.samples = [];
data.policyMean = [3];
data.policyCov = [3];
data.prefFeedback = [];
data.absFeedback = [];
data.failedExperiments = [];

if exist([data.fileName, '.mat'], 'file')
    error([data.fileName, '.mat already exists!'])
else
    save([data.fileName, '.mat'], 'data');
end
    
