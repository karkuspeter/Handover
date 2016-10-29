function savem = frame_reps_human()

addpath('../')
addpath('../gpml')

startup

N = 5;
episodes = 50;
epsilon = .5;
eta = 1;

mn = [-20.8452  -5.6532   -1.1623   -0.8214];
cn = diag([10, 10, 5, 5].^2);
dim = length(mn);

savem = [mn(:)', diag(cn)'.^.5, eta, 0, 0, 0, 0];

options = optimset('Algorithm','active-set');
options = optimset(options, 'GradObj','on');
options = optimset(options, 'Display', 'off');
S = [];
r = [];
logHyp = [];

for e = 1:episodes
    
    Smodel = S;
    rmodel = r;
    % sample
    for i = 1:N
        snew = mvnrnd(mn, cn);
        S = [S; snew];
        
        [res, t] = sim_invpendulum(snew(:)');
        f = figure(1); clf, subplot(3,1, [1, 2]);
        plot(t, res(:, 1:2:end));
        title(['Current gain: ', num2str(mn(:)')])
        subplot(3,1,3);
        if ~isempty(logHyp)
            [predRew, stdPredRew] = predictWithFullGPModel(logHyp, Smodel, rmodel, snew);
            xx = 0:.01:1;
            yy = exp(-.5*(predRew-xx).^2/stdPredRew^2)*.1;
            plot(xx, yy)
            ylabel('reward prediction')
        else
            predRew = NaN;
        end
        
        
        axis([0 1 0 .1]);
        [rew, yrew] = ginput(1);
        disp(['You gave ',num2str(rew), ' reward.'])     
%         rew = input(['What reward do you give? (prediction = ', num2str(predRew), '): ']);
        r = [r; rew];
    end
    
    
    % update reward model
    disp('Updating GP model...')
    logHyp = getFullGPModel(S, r, 200);
    disp('...ready!')
    
    rfull = [];
    Sfull = [];
    % generating new rewards on demand
    extraSamples = input('How many extra samples hould I generate: ');
    if extraSamples > 0
        Shat = mvnrnd(mn, cn, extraSamples);
        [rhat, stdrhat] = predictWithFullGPModel(logHyp, S, r, Shat);
        noExtra = 20;
        rewSamples = bsxfun(@plus, rhat , bsxfun(@times, randn(extraSamples, noExtra), stdrhat));
            
        rfull = reshape(rewSamples', [], 1);
        Sfull = reshape(repmat(Shat', noExtra, 1), dim, [])';
    else
        rfull = r;
        Sfull = S;
    end

    q = ones(length(rfull), 1)/length(rfull);
    
    % optimize
    objfun = @(eta) reps_dual(eta, rfull, q, epsilon);
    eta = fmincon(objfun, eta, -1, -.01, [], [], [], [], [], options);
    
    % update policy
    p = exp(rfull/eta)/sum(exp(rfull/eta));
    Dkl = sum(p.*log(p./q));
    iDkl = sum(q.*log(q./p));
    mn = Sfull'*p;
    cn = (Sfull-repmat(mn', size(Sfull, 1), 1))'*((Sfull-repmat(mn', size(Sfull, 1), 1)).*repmat(p, 1, dim));
    
    rexp = rfull(:)'*p;
    rvar = (rfull(:)-rexp)'*((rfull(:)-rexp).*p);
    
    % save data
    savem = [savem; mn(:)', diag(cn)'.^.5, eta, rexp, Dkl, iDkl, rvar];
end

% for i = 1:(size(savem, 1)-1)
% 	th = .5*atan(2*savem(i+1, 5) * savem(i+1, 3) * savem(i+1, 4)/(savem(i+1,3)^2 - savem(i+1, 4)^2));
% 	[ex, ey] = calculateEllipse(savem(i+1, 1), savem(i+1, 2), savem(i+1, 3), savem(i+1, 4),th *180/pi);
% 	plot([savem(i, 1); savem(i+1, 1)], [savem(i, 2), savem(i+1, 2)], 'w', 'Linewidth', 2), hold on;
% 	plot(ex, ey, 'w--')
% end

figure()
subplot(4,2,1), plot(savem(1:end, 1:2)), ylabel('mean')
subplot(4,2,2), plot(savem(2:end, 7)), ylabel('E[R]')
subplot(4,2,3), plot(savem(1:end, 3:4)), ylabel('std')
subplot(4,2,4), plot(savem(2:end, 6)), ylabel('eta')
subplot(4,2,5), plot(savem(2:end, 8)), ylabel('Dkl')
subplot(4,2,6), plot(savem(1:end, 9)), ylabel('iDkl')
subplot(4,2,7), plot(savem(2:end, 10).^.5), ylabel('Var[R]^{0.5}')
end

