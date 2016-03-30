function [rewards, stateActions] = simulateCartPole(x0, numTrajs, dt, tend, forceStd, K)

if nargin < 2
    K =  [ -1.7321
        -4.7158
        -101.9341
        -29.2456];
    dt = .01;
    tend = 2;
    forceStd = .4;
    numTrajs = 1;
elseif nargin < 3
     K =  [ -1.7321
        -4.7158
        -101.9341
        -29.2456];
    dt = .01;
    tend = 2;
    forceStd = .4;
elseif nargin < 6
    K =  [ -1.7321
        -4.7158
        -101.9341
        -29.2456];
end



Q = diag([3 .1 4 .4]); R = 1;
L = 1;
m = 1;
M = 3;

rewFunc = @(x,u) -sum((x * Q).* x, 2) - u.^2 * R;



% [A, B] = linearizeDynamics(objfun, zeros(4, 1), 0);
% K = lqr(A, B, Q, R);
% K = K(:);

time = [0:dt:tend];


for trajs = 1:numTrajs
    x = x0;
    for i = 1:length(time)
        
        u = -x(:)'*K(:) + randn(1) * forceStd;
        
        stateActions(i, :, trajs) = [x(:)', u];
        rewards(i, trajs) = rewFunc(x(:)', u);
        
        objfun = @(t, xcurr) cartpoleDiffEq(t, xcurr, L, m, M, u);
        [~, y] = ode45(objfun, [0 dt], x);
        
        x = y(end, :)';
        
    end
end
end


function dx = cartpoleDiffEq(t, x, L, m, M, u)
% x = [x, dx, th, dth]
% L length to mass
% m pole mass
% M cart mass
% u force in positive x direction
dx = zeros(4, 1);

dx(1) = x(2);
dx(2) = (u + m * sin(x(3))* (L * x(4)^2 - 9.81 * cos(x(3))))/...
    (M + m * sin(x(3))^2);
dx(3) = x(4);
dx(4) = (- u * cos(x(3)) - m*L*x(4)^2 * sin(x(3)) * cos(x(3)) + (M +m) * 9.81 * sin(x(3))) / ...
    (L*(M + m * sin(x(3))^2));
end

function dx = cartpoleDiffEqWithControlGain(t, x, L, m, M, K)
% x = [x, dx, th, dth]
% L length to mass
% m pole mass
% M cart mass
% u force in positive x direction
dx = zeros(4, 1);
u = -x'*K(:);
dx(1) = x(2);
dx(2) = (u + m * sin(x(3))* (L * x(4)^2 - 9.81 * cos(x(3))))/...
    (M + m * sin(x(3))^2);
dx(3) = x(4);
dx(4) = (- u * cos(x(3)) - m*L*x(4)^2 * sin(x(3)) * cos(x(3)) + (M +m) * 9.81 * sin(x(3))) / ...
    (L*(M + m * sin(x(3))^2));
end