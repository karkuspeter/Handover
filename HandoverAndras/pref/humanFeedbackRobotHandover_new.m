function [fPrior, prefs, context, samples] = humanFeedbackRobotHandover_new()
 

datasheet = csvread('Robot preference data - Sheet1.csv');

for i = 1:size(datasheet, 1)
   
    filename = ['dataRecord', num2str(datasheet(i, 1)), '.txt'];
    fPrior(i) = datasheet(i, 2);
    context(i) = datasheet(i, 3);
    samples(i, :) = datasheet(i, 4:10);
    prefOver(i, :) = datasheet(i, 11:end);
    
    % trajectory data
    jani = dlmread(filename);
    time = jani(:, 1) + jani(:, 2)/1e9 - jani(1, 1);
    forces = [smooth(jani(:, [3 ])), smooth(jani(:, [ 5 ])), smooth(jani(:, [ 7]))];
    torques =  [smooth(jani(:, [9 ])), smooth(jani(:, [ 11 ])), smooth(jani(:, [ 13]))];;
    pos = jani(:, [15:17]);
    absForce = sqrt(sum(forces.^2, 2));
    absTorque = sqrt(sum(torques.^2, 2));
%     
%     clf, subplot(3,1,1),
%     plot(time, (forces)), legend('x','y','z')
%     ylabel('forces')
%     hold on, plot(time, absForce, 'k')
%     subplot(3,1,2),
%     plot(time, (torques)), legend('x','y','z')
%     hold on, plot(time,absTorque, 'k')
%     ylabel('torques')
%     subplot(3,2,5),
%     plot(diff(absForce)./.02)
%     ylabel('jerk')
%     subplot(3,2,6),
%     plot(diff(absTorque)./.02)
%     ylabel('torque jerk')
    
    ixstart = find(absForce > 5); ixstart = ixstart(1);
    [fmax, ixmax] = max(absForce);
    ixend = find(torques(ixmax:end, 1) < 0); ixend = ixmax+ ixend(1);
    thandover = time(ixend)-time(ixstart);
    maxjerk = max(abs(diff(absForce)./.02));
    maxtorquejerk = max(abs(diff(absTorque)./.02));
    disp(['start/end/total time /reward: ', num2str(time(ixstart)), ' / ', num2str(time(ixend)), ' / ', num2str(thandover), ' /', num2str(fPrior(i))])
    disp(['MAX force/jerk: ', num2str(fmax), ' / ', num2str(maxjerk), ' / ', num2str(context(i))]);
    
    
    
end

% prefs
prefs = [];