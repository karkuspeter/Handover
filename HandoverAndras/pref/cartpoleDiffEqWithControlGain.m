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