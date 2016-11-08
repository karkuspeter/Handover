clear all, close all

addpath('../gpml')
addpath('../SPGP_dist')
startup % of GPML
addpath('/home/kupcsik/Dropbox/workonline/progs/matlab/svn_projects/REPS/matlab/');
addpath('/home/kupcsik/Dropbox/workonline/progs/matlab/svn_projects/REPS/matlab/MotorPrimitives');
addpath('/home/kupcsik/Dropbox/workonline/progs/matlab/svn_projects/REPS/matlab/Environments/4LinkBalancing');
addpath('/home/kupcsik/Dropbox/workonline/progs/matlab/svn_projects/REPS/matlab/evaluation/QuadLinkThrowState');
addpath('/home/kupcsik/Dropbox/workonline/progs/matlab/svn_projects/REPS/matlab/REPS');
addpath('/home/kupcsik/Dropbox/workonline/progs/matlab/svn_projects/REPS/matlab/evaluation');
properties.dt = 0.01;
properties.Ts = 2;
properties.Noise_std = 0;
properties.newSamples = 10;
properties.maxSamples = 20 * properties.newSamples;
properties.numRBF = 10;
environment = QuadPendulumThrowing(properties);
properties.useWeights = true;
properties.useGoal = false;
properties.useGoalVelocity = false;
properties.useTau = false;
properties.alphaZ = 3;    
properties.useInitPosition = false;
properties.useInitVelocity = false;
properties.tunePDgains = false;
properties.tuneReleaseTime = true;
properties.tunePIDGains = false;
properties.usePD = true;
properties.usePID = false;
properties.useLQ = false;
environment.fnCost = @task_throwing;
properties.dimState = 2;
properties.dimAction = 1 + 4 * properties.numRBF;
[ model, properties, trial ]    = initREPS(properties);
properties = createDMPGenericRewardPDtuned(properties, environment);
properties.drawFunctions = {@QuadLinkDrawFunction};
properties.dimState = 2;
properties.numIterations = 100;
properties.epsilonAction = .3;
properties.initSigma = .1;
properties.optionInitFunc = @QuadLinkThrowInit2stateWBeta;
[model, properties] = properties.initModel(environment, properties);  
model.BetaA = [
      6.7288  -31.4628
   -5.9392   23.2078
   21.8163   -7.8109
    2.4760  -23.0373
   -3.2004  -76.8050
    0.1083  -21.4141
  -25.7622   77.0693
   -1.7281  -95.9756
   -4.6100  -57.6641
   13.2376   51.8267
   -4.2466  -24.4433
   -0.2640   -8.9374
  -24.3546   17.5606
   -8.5055   32.1915
   17.7074  -85.9945
    1.0041  -80.9889
    5.8449  -30.7944
   20.5506   15.7873
   -3.8211  -72.6864
    5.7918  -13.9218
    2.0776   67.9060
    3.4221  -23.7532
   -4.8061  -34.5568
   12.0493   31.9732
   -2.0287   80.0993
   -9.5557   17.7829
   -1.5057  -56.4169
    0.8181 -126.3499
   16.0739  -38.9978
   20.3626  -33.4608
    1.1580    2.6433
   -3.5569   12.7185
  -22.6850   -4.7035
   10.3779  -68.0949
    6.9797   45.6995
  -22.5956   29.2454
   -1.9135    8.7847
    7.2765  -64.4205
  -10.6717  -10.6219
   14.3421   11.2635
    0.0050   -0.0273
];

% Extra parameters from simulation
PGains = properties.environment.PGains;
DGains = properties.environment.DGains;

% Prediction parameters
K = zeros(4, 8);
pg = -diag(PGains);
dg = -diag(DGains);
K(1, 1:2) = [pg(1), dg(1)];
K(2, 3:4) = [pg(2), dg(2)];
K(3, 5:6) = [pg(3), dg(3)];
K(4, 7:8) = [pg(4), dg(4)];

control_input = @(x, k, r) k*(x(:)-r(:));
iscov = @(S, k) k*S;
iicov = @(S, k) k*S*k';

