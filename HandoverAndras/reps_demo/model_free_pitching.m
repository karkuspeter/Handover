clear all, close all

% noisy pitch function
% a = -20 - 10*sqrt(3)
% b = 20 + 20*sqrt(3)
pitch = @(t)(-37.3205)* t.^2 + (54.6410)* t + 20 +...
    normrnd(0, 1, size(t));
% reward funtion
R = @(t) - (40 - pitch(t)).^2;

% plot in log scale
%figure;
%t = 0 : 0.005 : 2;
%plot(t,-log(-R(t)));

% for drawing
track_mu = zeros(21,1);
track_sigma = zeros(21,1);
% get data from experiment
N_real = 20;
mu = 1;
sigma = 0.1;
zeros_real = zeros(N_real, 1);
t_real = zeros_real + normrnd(mu, sigma, size(zeros_real));
r_real = R(t_real);
track_mu(1,1) = track_mu(1,1) + mu;
track_sigma(1,1) = track_sigma(1,1) + sigma;


for episode = 2:21
    % polynomial regression
    p = polyfit(t_real, r_real, 4);
    %r_est = polyval(p, 0:0.005:2);
    %plot(t_real, r_real, 'o', 0:0.005:2, r_est,'+');
    
    % draw artifical samples
    N_samples = 500;
    zeros_samples = zeros(N_samples, 1);
    t_samples = zeros_samples + rand(size(zeros_samples)) * 2;
    r_samples = polyval(p, t_samples);
    
    % compute reps dual
    options = optimset('Algorithm','active-set');
    options = optimset(options, 'GradObj','on');
    options = optimset(options, 'Display', 'off');
    eta = abs(mean(r_samples));
    epsilon = 0.5;
    objfun = @(eta) reps_dual(eta, r_samples, epsilon);
    eta = fmincon(objfun, eta, -1, -.01, [], [], [], [], [], options);
    
    p = exp(r_samples/eta) / sum(exp(r_samples/eta));
    mu = t_samples' * p;
    sigma = sqrt(sum(p .* (t_samples - mu) .*  (t_samples - mu)));
    N_new = 2;
    zeros_new = zeros(N_new, 1);
    t_new = zeros_new + normrnd(mu, sigma, size(zeros_new));
    t_real = [t_real', t_new']';
    r_real = [r_real', R(t_new)']';
    track_mu(episode,1) = track_mu(episode,1) + mu;
    track_sigma(episode,1) = track_sigma(episode,1) + sigma;
end

plot(1:21, track_mu, '-', 1:21, track_sigma, '+')

