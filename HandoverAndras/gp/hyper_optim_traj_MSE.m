function [f, df] = hyper_optim_traj_MSE(loghyp, xtrain, ytrain, traj, polpar)
% Calculates the gradient of the log-hyperparameters to minimize the mean squared error of trajectory prediction.
%
% Input:
% Output:

e = size(ytrain, 2);
d = size(xtrain, 2);

if and(or(size(loghyp, 2) == 1, size(loghyp, 1) == 1), e > 1)
    loghyp = reshape(loghyp, d+2, e);
end

hyp = exp(loghyp);
% hyp = loghyp;
% Init Mean Predict
for i = 1:e
    sigf = hyp(end-1,i); sign = hyp(end, i); w = hyp(1:end-2, i);
    K(:, :, i) = sigf^2 * exp(-.5*maha(xtrain, xtrain, diag(w.^-2)));
    L = chol(K(:, :, i)+eye(size(K, 1))*sign^2);
    mf(:, :, i) = eye(size(K, 1))/L/L'*ytrain(:, i);
end

figure(10), clf

for t = 1:size(traj, 3)
    trajPredLoc = traj(1, :, t);
    
    for i = 1:size(traj, 1)-1
        x = trajPredLoc(i, 1:e);
        u = -x*polpar';
        xquery = [x, u];
        
        meanpred = zeros(1, e);
        for j = 1:e
            w = hyp(1:end-2, j);
            sigf = hyp(end-1, j);
            k = sigf^2*exp(-.5*maha(xtrain, xquery, diag(w.^-2)));
            
            meanpred(j) = k'*mf(:, :, j);
        end
        
        trajPredLoc(i+1, 1:e) = trajPredLoc(i, 1:e) + meanpred;
        trajPredLoc(i+1, e+1:end) = trajPredLoc(i+1, 1:e)*polpar';
        
    end
    
    trajPred(:, :, t) = trajPredLoc;
    
    hold on, plot(trajPred(:, :, t))
end



N = size(traj, 3);
T = size(traj, 1);

f = 0;

for i = 1:size(traj, 3);

    f = f + sum(sum((traj(:, 1:e, i) - trajPred(:, 1:e, i)).^2));
    
end

grad = zeros(e*(d+2), e, T, N);
gradWeighted = zeros(e*(d+2), e, T, N);
gradW = zeros(e, T, N);

if nargout > 1
    
    % Precalculations
    for j = 1:e
        sigf = hyp(end-1, j);
        sign = hyp(end, j);
        w = hyp(1:end-2, j);
        W(:, :, j) = diag(w.^-2);
        
        K(:, :, j) = sigf^2 * exp(-.5*maha(xtrain, xtrain, W(:, :, j)));
        L = chol(K(:, :, j)+eye(size(K(:, :, j)))*sign^2); % Inverse: eye()/L/L';
        Kinv(:, :, j) = eye(size(L))/L/L';
        KinvY(:, :, j) = Kinv(:, :, j)*ytrain(:, j);
%         KinvY = L\L'\ytrain(:, j);
    end
    
    for i = 1:N % trajectory
        trajLoc = traj(:, :, i);
        predTrajLoc = trajPred(:, :, i);
        
        for t = 2:T % timestep 
            queryLoc = predTrajLoc(t-1, :);
            predTargetLoc = predTrajLoc(t, 1:e);
            realTargetLoc = trajLoc(t, 1:e);
            
            linDiff = bsxfun(@minus, xtrain, queryLoc);
            
            dMu_dTheta = grad(:, :, t-1, i);

            dQuery_dTheta = [dMu_dTheta'; (-polpar*dMu_dTheta')]';
            
            for j = 1:e % for each output dimension
                sigf = hyp(end-1, j);
                sign = hyp(end, j);
                w = hyp(1:end-2, j);
                
                beta = (W(:, :, j)*linDiff');
                
                k = sigf^2 * exp(-.5* maha(xtrain, queryLoc, W(:, :, j)));
                
                for jj = 1:e % for each GP model hyper-parameter vector
                    ixLoc = (jj-1)*(d+2)+1:jj*(d+2);
                    
                    if jj == j % calculate the j-dimension terms gradient w.r.t. j dimension hyperparameters
                        grad_sigf = -1/sigf * k' .* (dQuery_dTheta(ixLoc(end-1), :) * beta + (beta' * dQuery_dTheta(ixLoc(end-1), :)')');
                        Grad_sigf = 2/sigf * K(:, :, j);
                    else % calculate the j-dimension terms gradient w.r.t. NOT j dimension hyperparameters
                        grad_sigf = -.5 * k' .* (dQuery_dTheta(ixLoc(end-1), :) * beta + (beta' * dQuery_dTheta(ixLoc(end-1), :)')');
                    end
                    
                    for dd = 1:d % for each weight dimension
                        Wloc = zeros(d, d); Wloc(dd, dd) = w(dd)^-3;
                        if jj == j % calculate the j-dimension terms gradient w.r.t. j dimension hyperparameters
                            grad_w(dd, :) = .5 * k'.* ( -dQuery_dTheta(ixLoc(dd), :) * beta - (beta' * dQuery_dTheta(ixLoc(dd), :)')' -2 * maha(xtrain, queryLoc, Wloc)');
                            Grad_w(:, :, dd) = .5 * K(:, :, j) .* -2*maha(xtrain, xtrain, Wloc);
                        else % calculate the j-dimension terms gradient w.r.t. NOT j dimension hyperparameters
                            grad_w(dd, :) = .5 * k'.* ( -dQuery_dTheta(ixLoc(dd), :) * beta - (beta' * dQuery_dTheta(ixLoc(dd), :)')');                            
                        end
                    end
                    
                    grad_sign = zeros(1, length(k));
                    Grad_sign = eye(length(k)) * 2 * sign;
                    
                    grad_k = [grad_w; grad_sigf; grad_sign];
                    
                    grad(ixLoc, j, t, i) = dMu_dTheta(ixLoc, j) + grad_k * KinvY(:, :, j);
                    if jj == j % the gradK_gradTheta term is only nonzero for the actual output dimension
                        for dd = 1:d
                            grad(ixLoc(dd), j, t, i) = grad(ixLoc(dd), j, t, i) - k'*Kinv(:, :, j) * Grad_w(:, :, dd) * KinvY(:, :, j);
                        end
                        grad(ixLoc(end-1), j, t, i) = grad(ixLoc(end-1), j, t, i) - k'*Kinv(:, :, j) * Grad_sigf * KinvY(:, :, j);
                        grad(ixLoc(end), j, t, i) = grad(ixLoc(end), j, t, i) - k'*Kinv(:, :, j) * Grad_sign * KinvY(:, :, j);
                    end
                end
                % Finished with the first dimension of grad
                
            end
            % Finished with the second dimension of grad
            gradW(:, t, i) = predTargetLoc - realTargetLoc;
            gradWeighted(:, :, t, i) = bsxfun(@times, grad(:, :, t, i), gradW(:, t, i)');
        end
        % Finished with the third dimension of grad
        
        
    end
    % Finished with the fourth dimension of grad
    % turn to log derivate
    df = sum(sum(sum(gradWeighted, 4), 3), 2) .* reshape(hyp, [], 1);
end

% f, [df, loghyp]
% pause(.5)