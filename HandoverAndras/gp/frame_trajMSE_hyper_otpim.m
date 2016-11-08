clear all, close all
addpath('../gpml')
startup
% Pars
PD = [2.2532, 4.2274];
dt = .05;
tend = 2;

figure
xtrain = [];
ytrain = [];

for i = 1:25
    state = 2*(rand(1, 2)-.5);
    modelSeed = i;

    sim('noisyPendulum')

    % Init
    x = [simout.signal1.Data,  simout.signal2.Data];
    u =  simout.signal3.Data;
    t =  simout.signal4.Data;

    traj(:, :, i) = [x, u];

    hold on, plot(t, [x, u]);

    xtrain = [xtrain; [x(1:end-1, :), u(1:end-1)]];
    ytrain = [ytrain; diff(x)];

end

hyp_init = repmat(std(xtrain)'/10, 1, 2);
hyp_init = [hyp_init; std(ytrain); std(ytrain)/10];

logHypInit = getFullGPModel(xtrain, ytrain, 200, log(hyp_init));

startState = [.6, -.5];
figure,
for i = 1:25 
    state = startState;
    modelSeed = i;

    sim('noisyPendulum')

    % Init
    x = [simout.signal1.Data,  simout.signal2.Data];
    u =  simout.signal3.Data;
    t =  simout.signal4.Data;
    
    trajTest(:, :, i) = [x, u];

    hold on, plot(t, [x, u]);
end

figure,
trajTestMean = mean(trajTest, 3);
trajTestStd = std(trajTest, 0, 3);
plot(t, trajTestMean, 'LineWidth', 2), hold on
plot(t, trajTestMean + trajTestStd*2, '--', 'LineWidth', 2)
plot(t, trajTestMean - trajTestStd*2, '--','LineWidth', 2)

initStates = squeeze(traj(1, :, :))';

hyp = exp(logHypInit(:, 1));
W = diag(trajTestStd(2, 1).^-2);
simi = exp(-.5*maha(initStates(:, 1), trajTestMean(1, 1), W));



% [f, df] = hyper_optim_traj_MSE(reshape(logHypInit, [], 1), xtrain, ytrain, traj, PD)
% [f, df] = hyper_optim_traj_MSE_numerical(reshape(logHypInit, [], 1), xtrain, ytrain, traj, PD)

% [f] = hyper_optim_traj_MSE_numerical(log(hyp_init), xtrain, ytrain, traj, PD);
% disp(['Objective with initialization: ', num2str(f)])
% pause
% [f] = hyper_optim_traj_MSE_numerical(logHypInit, xtrain, ytrain, traj, PD);
% disp(['Objective with standard GP training: ', num2str(f)])
% % reshape(dfn, 5, 2)
% pause
% 
% 
% % hyp_opt = minimize(reshape(log(hyp_init), [], 1), @hyper_optim_traj_MSE, -200, xtrain, ytrain, traj, PD);
% hyp_opt = minimize(reshape(logHypInit, [], 1), @hyper_optim_traj_MSE_numerical, -100, xtrain, ytrain, traj, PD);