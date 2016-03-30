function [fPrior, prefs, context, samples, forceJerk, ixSlow] = humanFeedbackRobotHandover()

fPrior = [1 7
          4 8
          7 7
          10 8
          12 9
          15 5
          22 9
          23 5
          25 7
          26 4
          27 9
          28 7
          29 5
          32 7
          33 9
          36 3
          37 9 
          38 10
          40 4
          43 9
          44 10
          49 4
          51 7
          54 3
          55 9
          56 8
          61 4
          62 5
          64 6
          66 4
          67 6
          73 6
          76 5
          79 4
          80 3
          82 8
          83 7
          87 9
          92 4
          94 6
          97 3
          99 4];
    
      ixGood = ones(42, 1);
      ixGood(1:2) = 0;
      ixGood(10:15) = 0;
      fPrior(1:10, 1) = fPrior(1:10, 1) - 5;
      fPrior(11:end, 1) = fPrior(11:end, 1) - 5-10;
      fPrior = fPrior(ixGood == 1, :);
      
      % batch 8
      
      fPriorB8 = [101 4
          102 8
          104 7
          107 9
          108 8
          109 9
          118 7
          121 8
          122 5
          124 7
          127 10
          128 2
          129 9
          131 5
          132 6
          139 4
          141 8
          144 3
          149 4
          150 5];
      fPriorB8(:, 1) = fPriorB8(:, 1) - 15;
      fPriorB8(3:end, 1) = fPriorB8(3:end, 1) -1;
      
    
      fPrior = [fPrior; fPriorB8];
      
      fPriorB11 = [151 7;
          153 8
          156 8];
      fPriorB11(:, 1) = fPriorB11(:, 1) -16;
      
      fPriorB12 = [161 5
          165 5
          167 5];
      fPriorB12(:, 1) =  fPriorB12(:, 1) -16;
      
      fPriorB13 = [172 9
          174 8
          175 5
          176 6
          178 6
          181 2
          182 3
          ];
      fPriorB13(:, 1) =  fPriorB13(:, 1) -16;
      
      fPriorB14 = [183 6
          184 9
          186 4
          189 4
          190 5
          193 7
          196 8];
      fPriorB14(:, 1) =  fPriorB14(:, 1) -16;
      
          
      
      fPrior = [fPrior; fPriorB11; fPriorB12; fPriorB13; fPriorB14];
      
%      prefs = [1 2
%          2 3
%          4 3
%          4 5
%          5 6
prefs = [         7 6
         7 8
         8 9
         10 9
         10 11
         12 11
         12 13
         13 14
         15 14
         16 15
         16 17
         18 17
         18 19
         19 20
         21 20
         22 21
         22 23 
         24 23
         25 24 
         37 36
         38 37
         38 39
         39 40
         41 40
         42 40
         43 42
         44 43
         44 45
         46 45
         47 46
         47 48
         50 49
         51 52
         53 52
         55 54
         55 56
         57 56
         58 57
         58 59
         60 59
         60 61
         63 62
         64 65
         67 66
         67 68
         69 68
         69 70
         70 71
         72 71
         72 73
         74 73
         75 74
         76 77
         77 78
         81 80
         82 81
         82 83
         83 84
         85 84
         86 85
         87 86
         87 88
         88 89
         90 89
         90 91
         93 92
         95 94
         95 96
         96 97
         98 97
         98 99
         100 99];
     
     prefs = prefs - 5;
     prefs(20:end, :) = prefs(20:end, :) - 10;
 
     prefsB8 = [102 101
         104 105 
         106 105
         109 108
         109 110 
         111 110
         112 110
         112 111
         112 113
         109 114
         115 114
         115 116
         114 116 
         109 116
         116 117 
         112 117
         118 117
         118 119
         120 119
         121 120
         121 122
         122 123
         124 123
         124 125
         126 125
         127 126
         127 128
         129 128
         129 130
         130 131
         132 131
         133 132
         134 133
         134 135
         136 135
         136 137
         135 138
         140 139
         141 140
         141 142
         143 142
         143 144
         145 144
         146 145
         146 147
         148 147
         148 149];
     prefsB8 = prefsB8 - 15;
     prefsB8(2:end, :) = prefsB8(2:end, :) - 1;     
         prefs = [prefs; prefsB8];
         
         prefsB11 = [152 151
             152 153
             152 154
             155 154
             156 155
             156 157
             157 158
             159 158
             159 160]-16;
         prefs = [prefs; prefsB11];
         
         prefsB12 = [162 161
             162 163
             162 164
             164 165
             165 166
             168 167
             168 170
             171 170]-16;
         prefs = [prefs; prefsB12];
         
         prefsB13 = [172 174 
             172 175
             177 176
             179 178]-16;
         prefs = [prefs; prefsB13];
         
         prefsB14 = [184 183
             185 183
             184 186
             185 188
             188 189
             191 190
             192 191
             194 193
             195 194
             196 195]-16;
         prefs = [prefs; prefsB14];
     
     batchNums = [ones(6, 1); 2*ones(10, 1); 3*ones(10, 1); 4*ones(11, 1); 5*ones(15, 1); 6*ones(25, 1); 7*ones(25, 1)];
     
     ixSlow = [26 36 40 41 76 93 96 97 101]-16;
         
         samples = [];
         
parLoc = dlmread('Batch1ParametersLog.txt');
numSamples(1) = size(parLoc, 1);
samples = [samples; parLoc];
parLoc = dlmread('Batch2ParametersLog.txt');
numSamples(2) = size(parLoc, 1);
samples = [samples; parLoc];
parLoc = dlmread('Batch3ParametersLog.txt');
numSamples(3) = size(parLoc, 1);
samples = [samples; parLoc];
parLoc = dlmread('Batch4ParametersLog.txt');
numSamples(4) = size(parLoc, 1);
samples = [samples; parLoc];
parLoc = dlmread('Batch5ParametersLog.txt');
numSamples(5) = size(parLoc, 1);
samples = [samples; parLoc];
parLoc = dlmread('Batch6ParametersLog.txt');
numSamples(6) = size(parLoc, 1);
samples = [samples; parLoc];
parLoc = dlmread('Batch7ParametersLog.txt');
numSamples(7) = size(parLoc, 1);
samples = [samples; parLoc];


parLoc = dlmread('Batch8ParametersLog.txt');



ixOk = and(batchNums ~= 1 , batchNums ~= 4);

context = samples(ixOk, 1); 
context = [context; parLoc(12:end, 1)];
samples = samples(ixOk, 2:end);
samples = [samples; parLoc(12:end, 2:end)];

parLoc = dlmread('Batch11ParametersLog.txt');

context = [context; parLoc(:, 1)];
samples = [samples; parLoc(:, 2:8)];


parLoc = dlmread('Batch12ParametersLog.txt');

context = [context; parLoc(:, 1)];
samples = [samples; parLoc(:, 2:8)];

parLoc = dlmread('Batch13ParametersLog.txt');

context = [context; parLoc(:, 1)];
samples = [samples; parLoc(:, 2:8)];

parLoc = dlmread('Batch14ParametersLog.txt');

context = [context; parLoc(:, 1)];
samples = [samples; parLoc(:, 2:8)];


forceJerk = [151 32.3935    4.9333  126.5004 .26 .4
    152 23.9779    4.5971  124.4582 .24 .6
    153 40.4486    4.3862  369.3063 .4 .76
     154.0000   40.8547    4.6917  308.9926 .2 .6
 155.0000   37.2969    4.7869  229.2171 .28 .5
156.0000   26.7757    5.6165  154.1704 .2 .6
 157.0000   60.7891    7.6294  348.8349 .28 .6
158.0000   52.6315    4.5023  352.8979 .34  .66
159.0000   28.7860    5.3258  161.6021 .2 .4
160.0000   57.1037    5.0923  245.8038 .3 .6
161.0000   54.2182    5.6237  310.9801 .32 .6
  162.0000   23.6677    4.6240  116.2848 .26 .6
  163.0000   74.1263    6.2102  371.3722 .28 .52
   164.0000   20.5583    4.6305  166.9293 .16 .4
   165.0000   43.2905    4.4640  183.0907 .32 .64
   166.0000   58.7671    5.0104  214.4925 .34 .8
   167.0000   38.0958    4.4051  186.6092 .3 1
   168.0000   29.7884    5.3130  174.8802 .4 .6
    170.0000   28.5009    4.4621  147.6901 .34 .46
    171.0000   28.1238    4.1281  143.3254 .3 .47
    172.0000   29.6490    4.5540  140.7590 .28 .4
   174.0000   30.8710    4.6809  212.7435 .2 .46
175.0000   32.1306    4.4972  167.4177 .2 .5
 176.0000   51.5932    4.5609  270.7199 .4 .64
177.0000   32.4847    4.5485  176.5179 .32 .48
178.0000   67.9960    5.3723  374.8361 .3 .4
179.0000   57.1070    5.4943  292.5328 .26 .44
180.0000   39.9843    5.6691  208.9479 .46  .6
181.0000   68.1877    4.5980  411.6678 .32 1
182.0000   38.7412    5.0988  341.5501 .24 1
183  40.2385    5.3480  238.5892 .26 .45
    184 33.9426    4.7929  176.4525 .38 .5
    185 31.9056    4.6294  170.0453 .24 .6
    186 58.4843    5.1651  248.9848 .3 2
       187 44.9416    4.9534  215.5727 .28 1.8
       188 35.8779    4.3337  226.3139 .2 .4
     189 54.8579    6.3989  367.7241 .24 2
     190 44.5524    5.5221  288.5683 .3 1.8
     191  33.1498    4.6398  199.3326 .28 .4
192 66.3680    5.8978  341.7552 .24 .5
193 30.5004    4.8035  301.8045 .16 .5
194 41.7484    4.2783  337.2006 .2 .5
195 34.2479    6.6162  224.1692 .24 .4
196 48.9463    4.9685  317.1863 .16 .48];

forceJerk(:, 1) = forceJerk(:, 1) - 16;