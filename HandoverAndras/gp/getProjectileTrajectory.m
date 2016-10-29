function [trajx, trajy] = getProjectileTrajectory(x, dt);
%
% Calculates the 2D trajectory of a projectile, given its starting position and velocity
% Gives back useful data for GP model learning on the trajectory.
%
% Input:
%	x: the initial position and velocity of the projectile [1, 4]
%	dt: the time resolution [1]
%
% Output:
% 	trajx: the position and velocity trajectory until reaching y(1) = 0
% 	trajy: the dynamics of position and velcoity until reaching y(1) = 0

pos = x(1:2);
vel = x(3:4);

g = 9.81;
d = vel(1)/g*(vel(2) + sqrt(vel(2)^2 + 2*g*pos(2))) + pos(1);

T = (d - pos(1))/ vel(1);
t = 0:dt:(T+dt);
x = pos(1) + t*vel(1);
y = pos(2) + vel(2)*t - .5* 9.81*t.^2;

vx = ones(length(x), 1)*vel(1);
vy = ones(length(y), 1)*vel(2) - g*t(:);

trajx = [x(1:end)' y(1:end)' vx(1:end) vy(1:end)] + randn(length(x), 4).*repmat([.1, .1, .3, .3], length(x), 1)*dt;
trajy = diff(trajx);
trajx = trajx(1:end-1, :);

