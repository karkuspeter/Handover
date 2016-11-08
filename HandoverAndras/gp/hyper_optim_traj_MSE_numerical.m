function [f, df] = hyper_optim_traj_MSE_numerical(loghyp, xtrain, ytrain, traj, polpar)

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



N = size(traj, 1);
T = size(traj, 3);

f = 0;

for i = 1:size(traj, 3);

    f = f + sum(sum((traj(:, 1:e, i) - trajPred(:, 1:e, i)).^2));
    
end

func = @(par) hyper_optim_traj_MSE(par, xtrain, ytrain, traj, polpar);
df = numericalGradient(func, reshape(loghyp, 10, 1));