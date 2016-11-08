clear all, close all

data.fileName = 'ToyCannon';
data.epsilon = 0.5;
data.initSamples = 21;
data.updateSamples = 5;
data.samples = [0.5751    0.7159
    0.4918    0.8606
    0.6840    1.0001
    0.3667    1.0606
    0.3755    1.1609
    0.3757    0.9507
    0.5439    1.0997
    0.3188    0.9331
    0.3785    1.0146
    0.7677    1.0377
    0.5888    0.7335
    0.4165    0.6605
    0.4410    0.9798
    0.6351    1.0321
    0.2287    1.1671
    0.3404    0.8004
    0.4748    1.0128
    0.7860    0.8542
    0.5909    0.9317
    0.5561    1.1272];
data.policyMean = [0.5, 1];
data.policyCov = {[.2^2 0; 0 .2^2]};
data.prefFeedback = [2     1
     2     3
     4     3
     4     5
     6     5
     6     7
     8     7
     8     9
     9    10
    11    10
    11    12
    13    12
    13    14
    15    14
    15    16
    16    17
    19    20   ];
data.absFeedback = [  1     7
     3     2
     6     9
     7     3
     9     6
    10     3
    14     4
    15     9
    18     4];
data.failedExperiments = [];
data.hyp = [];

if exist(['./', data.fileName, '.mat'], 'file')
    error([data.fileName, '.mat already exists!'])
else
    save(['./',data.fileName, '.mat'], 'data');
end
    
