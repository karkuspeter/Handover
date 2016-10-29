clear all, close all

data.fileName = 'SimpleLearning';
data.epsilon = 0.5;
data.initSamples = 20;
data.updateSamples = 5;
data.samples = [1.5041
    0.8971
    2.9881
    1.6669
    2.6093
    1.1137
    3.9570
    5.6747
    0.4165
    1.1613
    1.9337
    2.6667
    1.6752
    0.5364
    2.6928
    5.4583
    3.3426
    1.6066
    4.4464
];
data.policyMean = [3];
data.policyCov = [3];
data.prefFeedback = [1     2
     3     2
     4     3
     4     5
     5     6
     6     7
     7     8
     9     8
     9    10
    11    10
    11    12
    13    12
    14    13
    14    15
    15    16
    17    16
    18    17
    18    19
    ];
data.absFeedback = [9     1
     4     7
     7    13
     5    17];
data.failedExperiments = [];

if exist([data.fileName, '.mat'], 'file')
    error([data.fileName, '.mat already exists!'])
else
    save([data.fileName, '.mat'], 'data');
end
    
