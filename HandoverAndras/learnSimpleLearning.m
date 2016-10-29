clear all, close all
addpath('./pref')
addpath('../Matlab_Network/')
addpath('./reps_demo/')
addpath('./gp/')

load('SimpleLearning.mat')
rew = @(x) 2*( .5*exp(- (x - 2) * 1 * (x -2)) + exp(- x * 4 * x)-.1174) ;


    % Sample from policy
    sample = mvnrnd(data.policyMean(end), data.policyCov(end), 1);
    
    % show the reward 
    disp(['Reward = ', num2str(rew(sample) + mvnrnd(0, 1)*.1)])
    
    % get feedback
    failed = input('Failed experiment? [y/N]: ', 's');
    if isempty(failed)
        failed = 'n';
    end
    
    
    if strcmp(failed, 'y')
        % ignore sample
        disp('ignoring sample')
        data.failedExperiments = [data.failedExperiments size(data.samples,1)];
    else
        % update data
        data.samples = [data.samples; sample];
        preferred = input('Is the current sample preferred over previous? [y/n/i]: ', 's');
        N = size(data.samples,1);
        if strcmp(preferred, 'y')
            if ~isempty(data.prefFeedback)
                data.prefFeedback = [data.prefFeedback; [N N-1]];
            else
                data.prefFeedback = [N N-1];
            end
        elseif strcmp(preferred, 'n')
            if ~isempty(data.prefFeedback)
                data.prefFeedback = [data.prefFeedback; [N-1 N]];
            else
                data.prefFeedback = [N-1 N];
            end
        end
        
        absfeedback = 11;
        while absfeedback > 10
            absfeedback = input('How do you rate the current sample (1-10)? [-1 to ignore]: ');
            if isempty(absfeedback)
                absfeedback = -1;
            end
        end
        if absfeedback > 0
            if ~isempty(data.absFeedback)
                data.absFeedback = [ data.absFeedback; [absfeedback, N]];
            else
                data.absFeedback = [absfeedback, N];
            end
        end
        
    end

    N = size(data.samples, 1);
   
    if (N >= data.initSamples) && (mod(N-data.initSamples, data.updateSamples) == 0)
        % update policy
        disp('updating policy')
    end
        


save('SimpleLearning', 'data')