clear all, close all
%
% addpath('../gp')
% addpath('..')
% addpath('~/svnprojects/ClassSystem/Helpers/')
% addpath('../DERIVESTsuite/')
% addpath('../teach-grad-hess/')


muContext = [5];
covContext = diag([1.5]);
rangeContext = [2 10];

angleStd = 1;
numTrajs = 1;
prefsForSamples = 3;
consecutivePref = 0;
initialSamples = 200;
newSamples = 20;
epsilon = .5;
policyUpdates = 200;
ridge = 1e-4;
updateSamples = 500;
controlContexts = 30;
absFeedbackFrequency = 3; % evalutes every #... sample
siga = .1;
deterministicEvaluation = 1; % always take the most likely controller to check performance
hyperRestarts = 8;
numIterations = 20;
useCmaes = 1;
cmaesSteps = 150;

saveFileName = ['PrefToyCannonSimple_REPS', num2str(prefsForSamples), '_Consec', num2str(consecutivePref), '_AbsFreq', num2str(absFeedbackFrequency), '_InitNewSamples'...
    ,num2str(initialSamples), '_', num2str(newSamples), '_eps', num2str(round(100*epsilon)), '_siga', num2str(round(siga)), '.mat'];

% try
%     matlabpool 4
% catch
% end

rng('default'); rng(666);
% xcheck = mvnrnd(muContext, covContext, controlContexts);
xcheck = rand(controlContexts, 1) * (rangeContext(2)-rangeContext(1))+rangeContext(1);

maxUpdateSamples = 300;
newSampless = [20 40 60];
epsilons = [ .5 .75 1];

for ns = 1:length(newSampless)
    for es = 1:length(epsilons)
        
        newSamples = newSampless(ns);
        epsilon = epsilons(es);
        
        saveFileName = ['RepsToyCannon_NewSamp', num2str(newSamples), '_eps', num2str(round(100*epsilon)),'.mat'];

        for iter = 1:numIterations
            
            rng('default')
            rng(iter)
            
            muK =  [pi/4];
            covK = diag([.3; ].^2);
            
            samples = [];
            rewards = [];
            contexts = [];
            prefs = [];
            fPrior = [];
            
            %% Get initial samples
            for i = 1:initialSamples
                K = mvnrnd(muK, covK);
                %         x0 = mvnrnd(muContext, covContext);
                x0 = rand(1)*(rangeContext(2)-rangeContext(1))+rangeContext(1);
                %         [r, xu] = simulateCartPole([x0(1) 0 x0(2) 0]', 1, dt, tend, forceStd, K);
                r = simulateToyCannon_SimplePolicy(K, x0);
                
                samples = [samples; K(:)'];
                contexts = [contexts; x0(:)'];
                rewards = [rewards; r];
                
                
            end
            
            A = zeros(1, 1);
            a = muK;
            
            
            disp(' ')
            disp(' checking initial performance ...')
            disp(' ')
            
            % check initital performance
            for i = 1:controlContexts
                
                x0 = xcheck(i, :);
                if deterministicEvaluation
                    K = a + A*x0(:);
                else
                    K = mvnrnd(a + A*x0(:), covK);
                end
                r = simulateToyCannon_SimplePolicy(K, x0);
                
                
                rcheck(i) = r;
            end
            
            RMeanSave(iter, 1) = mean(rcheck);
            RStdSave(iter, 1) = std(rcheck);
            
            
            
            %% Main loop
            for iterPolicyUpdate = 1:policyUpdates
                
                disp(' ')
                disp('New iteration, sampling... ')
                disp(' ')
                % Sample
                for i = 1:newSamples
                    
                    x0 = mvnrnd(muContext, covContext);
                    K = mvnrnd(a + A*x0(:), covK);
                    
                    r = simulateToyCannon_SimplePolicy(K, x0);
                    
                    samples = [samples; K(:)'];
                    contexts = [contexts; x0(:)'];
                    rewards = [rewards; r];
                    
                end
                
                ixUpdate = max(1, length(rewards)-maxUpdateSamples+1):length(rewards);
                
                
                disp(' ')
                disp('Updating policy... ')
                disp(' ')
                % Update policy
                Phi = [ones(length(ixUpdate), 1), contexts(ixUpdate, :), contexts(ixUpdate, :).^2];
                repsobjfun = @(etatheta) reps_dual_state(etatheta, epsilon, rewards(ixUpdate), Phi);
                options = optimoptions('fminunc','GradObj','on','Algorithm','trust-region', 'Display', 'off');
                etatheta = [log(1); randn(3, 1)];
                [etatheta, ~, flag] = fminunc(repsobjfun, etatheta(:), options);%, diag([-1; ones(3, 1)]), [-.0001;Inf; Inf; Inf], [], [], [], [], [], options);
                disp(['fminunc flag: ', num2str(flag)])
                
                eta = exp(etatheta(1));
                theta = etatheta(2:end);
                
                V = Phi*theta;
                
                p = exp((rewards(ixUpdate) - V)/eta);
                p = p/sum(p);
                
                if ~isnan(sum(p))
                    
                    
                    dKL = sum(p.*log(p./(ones(length(p), 1)/length(p))));
                    if (dKL - epsilon)^2 > .01
                        warning('real dKL is far from epsilon')
                    end
                    
                    CC = [ones(length(ixUpdate), 1),  contexts(ixUpdate, :)];
                    aA = (CC' * diag(p) * CC)\CC'*diag(p)*samples(ixUpdate, :);
                    
                    a = aA(1, :)';
                    A = aA(2:end, :)';
                    
                    dummy = samples(ixUpdate, :) - bsxfun(@plus, a , A*contexts(ixUpdate, :)')';
                    covK = bsxfun(@times, p, dummy)'*dummy / sum(p);
                else
                    warning('NaN p, ignoring policy update')
                end
                
                disp(' ')
                disp('Checking learning performance ...')
                disp(' ')
                % check performance
                for i = 1:controlContexts
                    
                    x0 = xcheck(i, :);
                    if deterministicEvaluation
                        K = a + A*x0(:);
                    else
                        K = mvnrnd(a + A*x0(:), covK);
                    end
                    
                    r = simulateToyCannon_SimplePolicy(K, x0);
                    
                    rcheck(i) = r;
                end
                
                
                disp([' ============= UPDATE #', num2str(iterPolicyUpdate), ' ============'])
                disp(' ')
                disp([' Fixed states reward mean/std: ', num2str(mean(rcheck)), '/', num2str(std(rcheck)) ]);
                disp([' Initial Fixed states reward mean/std: ', num2str(RMeanSave(1)), '/', num2str(RStdSave(1)) ]);
                disp([' Exploration stds: ', num2str(diag(covK)'.^.5) ]);
                disp(' ')
                disp([' ====================================='])
                
                RMeanSave(iter, iterPolicyUpdate+1) = mean(rcheck);
                RStdSave(iter, iterPolicyUpdate+1) = std(rcheck);
                
                R.mean = RMeanSave;
                R.std = RStdSave;
                R.policyA(:, :, iter) = A;
                R.policya(:, iter) = a;
                
                
                
                
                %                 figure(1)
                %                 clf
                %                 subplot(2,1,1),
                %                 plot(RMeanSave)
                %                 subplot(2,1,2),
                %                 plot(RStdSave)
                %                 drawnow
                
            end
            save(saveFileName, 'R');
        end
    end
end
