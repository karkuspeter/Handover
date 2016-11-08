function frame_hyper_optim_KLdiv_Traj()
clear all, close all

% matlabpool % in case you have parallel computing toolbox

% Initial hyperparameters
loghyp = log([0.1, 0.1, 0.1, 1, .01; 0.1, 0.1, 0.1, 1, .01])';

% Control input noise
torque_noise = .2;
% Control gain
K = -[3, 1];
% dt
dt = .1;
% Simulation yime length
tend = 2;
% Number of observed trajectories
notrajs = 10;

% The differential equation
diffx_int2 = @(x, k) dt*[x(2), cos(x(:)'*x(:)) +  x(:)'*k(:) + randn(1)*torque_noise]';
diffx_loc = @(x) diffx_int2(x, K);

tm = [0:dt:tend];
T = length(tm);
xsave = [];
usave = [];
mnsave = zeros(length(tm), notrajs, 2);

trainx = [];
trainy = [];

figure,
% Simulation
for e = 1:notrajs
    x = mvnrnd([.5 -.5],diag([.1 .1]/5))'; % <-- Initial state is drawn
    y = zeros(length(tm), 2);
    y(1, :) = x;
    for t = 2:length(tm)
        x = x + diffx_loc(x);
        y(t, :) = x';
    end
    xsave = [xsave; y];
    usave = [usave; y*K'];
    
    % GP training data
    trainx = [trainx; y(1:end-1, :) y(1:end-1, :)*K']; 
    trainy = [trainy; y(2:end, :) - y(1:end-1, :)];
    
    % Save simulation data in 3D: (time point, trajectory, state dimension)
    mnsave(:, e, 1) = y(:, 1);
    mnsave(:, e, 2) = y(:, 2);
    
    % Plot results
    plot(tm, [y, y*K']), hold on
end

figure,
subplot(2,1,1), plot(tm, squeeze(mean(mnsave, 2))), ylabel('state means')
subplot(2,1,2), plot(tm, squeeze(mean(mnsave.^2, 2))-squeeze(mean(mnsave, 2)).^2), ylabel('state variances')

% Initial feature (Squared only for now)
phi_init = [mean(squeeze(mnsave(1, :, :)), 1)'; mean(squeeze(mnsave(1, :, :)), 1)'.^2];
% phi_init = mean(squeeze(mnsave(1, :, :)), 1)';
th = randn(T, length(phi_init))*.1; % <--- initial theta 
th_loghyp_init = [reshape(th', [], 1); reshape(loghyp, [], 1)]; % <-- initial loghyperparameter-theta composite

objfun = @(th_loghyp) hyper_optim_kldiv(th_loghyp, mnsave, trainx, trainy, T, notrajs); % gives back the dual function value and its gradient w.r.t. th_loghyp
objfun = @(th_loghyp) gradientFrame(th_loghyp, mnsave, trainx, trainy, T, notrajs); % gives back the dual function value and its gradient w.r.t. th_loghyp

% Compare numerical and analytic gradients, they should be the same
% numGrad = numericalGradient(objfun, th_loghyp_init);
% [~, analGrad] = objfun(th_loghyp_init);
% [numGrad, analGrad]

trainy = bsxfun(@times, randn(size(trainy)), [.001, .01]/2) + trainy;

th_loghyp_opt = getFullGPModelGPU(trainx, trainy, -200, 0);
th_loghyp_opt
% Perform optimization
th_loghyp_opt = minimize(th_loghyp_init, objfun, -200);
keyboard
end

function [f, df] = hyper_optim_kldiv(th_loghyp, mnsave, trainx, trainy, T, notrajs)
% This function computes the dual function and its gradient
% The current inputs are only for the time being, we should simplify it
% later!

% Some parameters, needs to be automated later
thetaFactor = .01;
xdim = 2;
udim = 1;
fulldim = xdim + udim;

% Reformating theta and loghyp from the parameter vector
th = reshape(th_loghyp(1:end-(xdim*(fulldim+2))), 2*xdim, [])';
% th = reshape(th_loghyp(1:end-(xdim*(fulldim+2))), xdim, [])';
loghyp = reshape(th_loghyp(end-(xdim*(fulldim+2))+1:end), fulldim+2, []);

% Initial feature (Squared for now) (should be passed to this function?)
phi_init = [mean(squeeze(mnsave(1, :, :)), 1)'; mean(squeeze(mnsave(1, :, :)), 1)'.^2];
% phi_init = [mean(squeeze(mnsave(1, :, :)), 1)'];

% Get the prediction means and variances at each training input, calculate
% their gradient too (grad_mean, grad_var) w.r.t. the log-hyperparameters
% Note: to avoid GP model overfitting, we leave always the query
% input-target pair from the training database. This way the query target
% remains hidden. I'm not sure if it's necessary, but it shouldn't hurt :)
[pred_mean, pred_var, grad_mean, grad_var] = getGradientsLeaveOneOut(loghyp, trainx, trainy, 1); 

% reshape the predictions such that we have 3 dimensions: (time, trajectory,
% state dimension)
pred_mean = reshape(pred_mean, [T-1, notrajs, xdim]); 
pred_mean = shiftdim(pred_mean, 1);

dummyMnsave = shiftdim(mnsave, 1);
err_mean = (dummyMnsave(:, :, 2:end) - pred_mean).^2;
err_mean = shiftdim(err_mean, 1);
figure(6),clf
for i = 1:notrajs
    plot(err_mean(:, :, i)'), hold on
end

pred_var = reshape(pred_var, [T-1, notrajs, xdim]);
pred_var = shiftdim(pred_var, 1);
% Reshape gradients according to the prediction mean and variance
grad_mean = reshape(grad_mean, [fulldim+2, T-1, notrajs, xdim]);
grad_var = reshape(grad_var, [fulldim+2, T-1, notrajs, xdim]);

dummySave = [];

% Now, compute the Z factors and the gradients of theta
for i = 1:T
    % Compute feature (squared for now)
    phi(:, :, i) = [squeeze(mnsave(i, :, 1:2))' ; squeeze(mnsave(i, :, 1:2))'.^2];
%     phi(:, :, i) = [squeeze(mnsave(i, :, 1:2))'];
    % Compute expected feature of the successor state (we don't need
    % E[phi(x_{T+1})] )
    if i < T
        expPhi(:, :, i) = [pred_mean(:, :, i)' ; pred_mean(:, :, i)'.^2 + pred_var(:, :, i)'];
%         expPhi(:, :, i) = [pred_mean(:, :, i)' ];
    end
    if i < T

        % To compute the gradients we don't directly get Z, but the
        % argument of Z, which we denote by ZZ. To avoid
        % \frac{sum(exp(...))}{\sum(exp(...)} divisions and their numerical
        % problems, we rather obtain the log-sum-exp trick:
        % instead of sum_i(exp(ZZ(i))*a/sum(exp(ZZ))) we use
        % exp( log(exp(ZZ(i))*a ) - log(sum(exp(ZZ))) ) = 
        % = exp( ZZ(i) + log(a) - log(sum(exp(ZZ))) ) = 
        % = exp( ZZ(i) + log(a) - (A + log(sum(exp(ZZ-A)))))
        % where A = max(ZZ(i))

        ZZ(i, :) = th(i, :) * phi(:, :, i) - th(i+1, :) * expPhi(:, :, i);

        if i == 1
            % This is what we want to compute
            % grad_th(i, :) = phi_init(:)' -  (phi(:, :, i)*Z(i,:)'/(sum(Z(i, :)) + 1e-5))'; 
            % This is how we do it
            A = max(ZZ(i, :));
            grad_th(i, :) = phi_init(:)' - ...
                sum(bsxfun(@times, exp(ZZ(i, :) - (A + log(sum(exp(ZZ(i, :) - A))))), phi(:, :, i)), 2)' ...
                + thetaFactor *th(i, :);
            
            dummySave(i, :, 1) = phi_init(:)';
            dummySave(i, :, 2) = sum(bsxfun(@times, exp(ZZ(i, :) - (A + log(sum(exp(ZZ(i, :) - A))))), phi(:, :, i)), 2)';
            
            dummySave(i, :, 3) = sum(bsxfun(@minus, phi_init(:) , phi(:, :, i)).^2, 2)';
        else
            % This is what we want to compute
            % grad_th(i, :) = (expPhi(:, :, i-1)*Z(i-1, :)'/(sum(Z(i-1, :)) + 1e-5))' - (phi(:, :, i)*Z(i, :)'/(sum(Z(i, :)) + 1e-5))'
            % This is how we do it
            Ai = max(ZZ(i, :));
            Aim1 = max(ZZ(i-1, :));
            grad_th(i, :) = sum(bsxfun(@times, exp(ZZ(i-1, :) - (Aim1 + log(sum(exp(ZZ(i-1, :) - Aim1))))), expPhi(:, :, i-1)), 2)' -  ...
                sum(bsxfun(@times, exp(ZZ(i, :) - (Ai + log(sum(exp(ZZ(i, :) - Ai))))), phi(:, :, i)), 2)' ...
                + thetaFactor *th(i, :);
            
            dummySave(i, :, 1) = sum(bsxfun(@times, exp(ZZ(i-1, :) - (Aim1 + log(sum(exp(ZZ(i-1, :) - Aim1))))), expPhi(:, :, i-1)), 2)';
            dummySave(i, :, 2) = sum(bsxfun(@times, exp(ZZ(i, :) - (Ai + log(sum(exp(ZZ(i, :) - Ai))))), phi(:, :, i)), 2)';
            
            
            dummySave(i, :, 3) = sum((expPhi(:, :, i-1) - phi(:, :, i)).^2, 2)';
        end
    else % i == T

        ZZ(i, :) = th(i, :) * phi(:, :, i);
        
        Ai = max(ZZ(i, :));
        Aim1 = max(ZZ(i-1, :));

        grad_th(i, :) = sum(bsxfun(@times, exp(ZZ(i-1, :) - (Aim1 + log(sum(exp(ZZ(i-1, :) - Aim1))))), expPhi(:, :, i-1)), 2)' -  ...
                sum(bsxfun(@times, exp(ZZ(i, :) - (Ai + log(sum(exp(ZZ(i, :) - Ai))))), phi(:, :, i)), 2)' ...
                + thetaFactor *th(i, :);
            
            dummySave(i, :, 1) = sum(bsxfun(@times, exp(ZZ(i-1, :) - (Aim1 + log(sum(exp(ZZ(i-1, :) - Aim1))))), expPhi(:, :, i-1)), 2)';
            dummySave(i, :, 2) = sum(bsxfun(@times, exp(ZZ(i, :) - (Ai + log(sum(exp(ZZ(i, :) - Ai))))), phi(:, :, i)), 2)';
            
            dummySave(i, :, 3) = sum((expPhi(:, :, i-1) - phi(:, :, i)).^2, 2)';
    end
end

if any(any(isnan(grad_th)))
    keyboard
end

figure(7), clf
plot(dummySave(:, :, 3),'-')
drawnow

figure(8), clf
plot(dummySave(:, :, 1), '--')
hold on, plot(dummySave(:, :, 2), '-')
drawnow

% pause
% ----------------
% Testing theta gradient 
% ----------------
% keyboard
% func = @(th1) -log(sum(exp(th1(:)'*phi(:, :, 1) - th(2, :)*expPhi(:, :, 1)))) + th1(:)'*phi_init;
% [numericalGradient(func, th(1, :)') grad_th(1, :)']
% func = @(th2) -log(sum(exp(th2(:)'*phi(:, :, 2) - th(3, :)*expPhi(:, :, 2)))) -log(sum(exp(th(1, :)*phi(:, :, 1) - th2(:)'*expPhi(:, :, 1))));
% [numericalGradient(func, th(2, :)') grad_th(2, :)']
% ----------------
% Testing GP model gradients
% ----------------
% checkDim = 1;
% [predMean, predVar, analMeanGrad, analVarGrad] = getGradientsLeaveOneOut(loghyp(:, checkDim), trainx, trainy(:, checkDim), 1);
% for i = 1:5
%     
%     funcMean = @(lhyp) dummyFuncMean(lhyp, trainx, trainy(:, checkDim), i);
%     funcVar = @(lhyp) dummyFuncVar(lhyp, trainx, trainy(:, checkDim), i);
%     numGradMean = numericalGradient(funcMean, loghyp(:, checkDim));
%     numGradVar = numericalGradient(funcVar, loghyp(:, checkDim));
%     
%     [[numGradMean; numGradVar], [analMeanGrad(:, i); analVarGrad(:, i)]]
% end
%    keyboard

% Computing the expected successor state feature gradient for each time
% step for each trajectory for each GP model
for j = 1:(T-1) % time steps
    for k = 1:notrajs % trajectories
        grad_expphi_loc = zeros(2*xdim, fulldim+2, xdim);
%         grad_expphi_loc = zeros(xdim, fulldim+2, xdim);
        for i = 1:xdim % GP models
            grad_mean_loc = squeeze(grad_mean(:, j, k, i)); 
            grad_var_loc = squeeze(grad_var(:, j, k, i));
            
            grad_expphi_loc(i, :, i) = grad_mean_loc'; % gradient of mean(x)
            grad_expphi_loc(i+xdim, :, i) = 2*pred_mean(k, i, j)*grad_mean_loc' + grad_var_loc'; % gradient of mean(x^2) = mean(x)^2 + var(x)
            
            grad_loghyp(k, :, i) = th(j+1, :)*grad_expphi_loc(:, :, i); % theta_t^T * grad(E[phi])/grad(w_i)
        end
    end
    
    for i = 1:xdim
        A = max(ZZ(j, :));
        grad_loghyp_final(j, :, i) = sum(bsxfun(@times, exp(ZZ(j, :) - (A + log(sum(exp(ZZ(j, :) - A))))), grad_loghyp(:, :, i)'), 2)';
%         grad_loghyp_final(j, :, i) = Z(j, :) * grad_loghyp(:, :, i)/(sum(Z(j, :))+1e-5);
    end
end
clear grad_loghyp
grad_loghyp = reshape(squeeze(sum(grad_loghyp_final, 1)), [], 1);

figure(9), clf 
% Plotting predictive mean and variance for all trajectories
for i = 1:notrajs
    subplot(2,1,1)
    plot(squeeze(pred_mean(i, :, :))'), hold on, ylabel('Pred mean')
    subplot(2,1,2)
    plot(squeeze(pred_var(i, :, :))'), hold on, ylabel('Pred var')
end

figure(10), clf
subplot(2,1,1), plot(grad_th); ylabel('gradient of thetas vs time')
subplot(2,1,2), plot(reshape(grad_loghyp, [], xdim)); ylabel('grad loghyp')

grad_th = reshape(grad_th', [], 1);

df = [grad_th; grad_loghyp];
f = th(1, :)*phi_init ;
% again, using the logsumexp trick for numerical accuracy
for i = 1:T
    A = max(ZZ(i, :));
    f = f - (A + log( sum(exp(ZZ(i, :) - A))))+ thetaFactor * th(i, :) * th(i, :)' ;
end

title(['Obj fun: ', num2str(f)])
drawnow

% disp('exp(ZZ)')
% exp(ZZ)

% 
% disp('Hyper-parameters')
% exp(loghyp)
disp('theta')
th
% % disp('Grad Hyperparameters')
% reshape(grad_loghyp, [], 2)
% disp('Grad theta')
% reshape(grad_th, 4, [])'
disp('Norm theta')
disp(norm(th))

% As we are using a minimization function, we change the signs
f = -f;
df = -df;

end

function [f, df] = gradientFrame(th_loghyp, mnsave, trainx, trainy, T, notrajs)
% to compare numerical and analytic gradients

[f, df] = hyper_optim_kldiv(th_loghyp, mnsave, trainx, trainy, T, notrajs); 

% objfun = @(par) hyper_optim_kldiv(par, mnsave, trainx, trainy, T, notrajs); % gives back the dual function value and its gradient w.r.t. th_loghyp
% numGrad = numericalGradient(objfun, th_loghyp);
% disp(['Norm(numgrad - analgrad): ', num2str(norm(numGrad-df))])

end

function [pred_mean] = dummyFuncMean(loghyp, trainx, trainy, i)
    [pred_mean, ~] = getGradientsLeaveOneOut(loghyp, trainx, trainy, 1);
    pred_mean = pred_mean(i);    
end

function [pred_var] = dummyFuncVar(loghyp, trainx, trainy, i)
    [~, pred_var] = getGradientsLeaveOneOut(loghyp, trainx, trainy, 1);
    pred_var = pred_var(i);
end

function [grad_mean, grad_var] = getGradients(loghyp, x, y)

e = size(y, 2);
[n, d] = size(x);

if and(or(size(loghyp, 2) == 1, size(loghyp, 1) == 1), e > 1)
    loghyp = reshape(loghyp(:), d+2, e);
end

hyp = exp(loghyp);

grad_mean = zeros(d+2, n, e);
grad_var = zeros(d+2, n, e);

for i = 1:e
    sigf = hyp(end-1,i); sign = hyp(end, i); w = hyp(1:end-2, i);
    K = sigf^2 * exp(-.5*maha(x, x, diag(w.^-2)));
    L = chol(K+eye(n)*sign^2);
    
    alpha = eye(n)/L/L'*y(:, i);
    
    for j = 1:d
        Wloc = zeros(d);
        Wloc(j, j) = w(j)^-3;
        gradLoc(:, :, j) = K.*maha(x, x, Wloc);
    end
    
    gard_mean_loc = zeros(d+2, n);
    gard_var_loc = zeros(d+2, n);
    
    parfor ii = 1:n
        
        k = sigf^2 * exp(-.5*maha(x, x(ii, :), diag(w.^-2)));
        beta = k'/L/L';
        
        % grad w
        for j = 1:d+2
            if j <= d
                Wloc = zeros(d);
                Wloc(j, j) = w(j)^-3;
                gradK = k.* maha(x, x(ii, :), Wloc);
                grad_mean_loc(j, ii) = (gradK' - beta*gradLoc(:, :, j)) * alpha;
                grad_var_loc(j, ii) = - 2*gradK' * beta' + beta*gradLoc(:, :, j)* beta';
            end
            
            if j == d+1
                % grad sigf
                grad_mean_loc(j, ii) = 2/sigf*(k' - beta*K)*alpha;
                grad_var_loc(j, ii) = 2*sigf - 2*2/sigf * k'*beta' + 2/sigf * beta * K * beta';
            end
            
            if j == d+2
                % grad sign
                grad_mean_loc(j, ii) = - 2*sign * beta * alpha;
                grad_var_loc(j, ii) =   2*sign * (beta * beta') + 2*sign;
            end
            
        end
    end
    
    grad_mean(:, :, i) = bsxfun(@times, grad_mean_loc, hyp(:, i));
    grad_var(:, :, i) = bsxfun(@times, grad_var_loc, hyp(:, i));
end



end

function [pred_mean, pred_var, grad_mean, grad_var] = getGradientsLeaveOneOut(loghyp, x, y, constPrior)
% This is an unoptimized implementation, might be slow, but we can put it
% on GPU and then it will be fast (hopefully :)). Nevertheless it works fine.

if nargin < 4 
    constPrior = 1; % add the query to the predictive mean
end

e = size(y, 2);
[n, d] = size(x);

if and(or(size(loghyp, 2) == 1, size(loghyp, 1) == 1), e > 1)
    loghyp = reshape(loghyp(:), d+2, e);
end

hyp = exp(loghyp);

pred_mean = zeros(n, e);
pred_var = zeros(n, e);

grad_mean = zeros(d+2, n, e);
grad_var = zeros(d+2, n, e);

for i = 1:e
    
    gard_mean_loc = zeros(d+2, n);
    gard_var_loc = zeros(d+2, n);
    parfor ii = 1:n
        
        ixLoc = [1:ii-1, ii+1:n];
%         ixLoc = [1:n];
    
        sigf = hyp(end-1,i); sign = hyp(end, i); w = hyp(1:end-2, i);
        K = sigf^2 * exp(-.5*maha(x(ixLoc, :), x(ixLoc, :), diag(w.^-2)));
        L = chol(K+eye(length(ixLoc))*sign^2);
        
        k = sigf^2 * exp(-.5*maha(x(ixLoc, :), x(ii, :), diag(w.^-2)));
        k = k(:);
        beta = k'/L/L';
    
        alpha = eye(length(ixLoc))/L/L'*y(ixLoc, i);
        
        pred_mean(ii, i) = k'*alpha;
        if constPrior 
            pred_mean(ii, i) = pred_mean(ii, i) + x(ii, i);
        end
        pred_var(ii, i) = sigf^2  - beta*k ; %+ sign^2;
    
        % grad w
        for j = 1:d+2
            if j <= d
                Wloc = zeros(d);
                Wloc(j, j) = w(j)^-3;
                gradK = k .* maha(x(ixLoc, :), x(ii, :), Wloc);
                gradLoc = K.*maha(x(ixLoc, :), x(ixLoc, :), Wloc);
                grad_mean_loc(j, ii) = (gradK' - beta*gradLoc) * alpha;
                grad_var_loc(j, ii) = - 2*gradK' * beta' + beta*gradLoc* beta';
            end
            
            if j == d+1
                % grad sigf
                grad_mean_loc(j, ii) = 2/sigf*(k' - beta*K)*alpha;
                grad_var_loc(j, ii) = 2*sigf - 2*2/sigf * k'*beta' + 2/sigf * beta * K * beta';
            end
            
            if j == d+2
                % grad sign
                grad_mean_loc(j, ii) = - 2*sign * beta * alpha;
                grad_var_loc(j, ii) =   2*sign * (beta * beta') ; %+ 2*sign;
            end
            
        end
    end
    
    grad_mean(:, :, i) = bsxfun(@times, grad_mean_loc, hyp(:, i));
    grad_var(:, :, i) = bsxfun(@times, grad_var_loc, hyp(:, i));
end

end

