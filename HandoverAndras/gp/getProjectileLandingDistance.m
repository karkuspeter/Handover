function d = getProjectileLandingDistance(x);
% 
% Calculates the landing distance (y = 0) of a projectile given the initial position and velocity (in 2D).
%
% Input:
%	x: pos and vel [1, 4]
%
% Output:
%	d: distnace [1]
pos = x(1:2);
vel = x(3:4);
g = 9.81;
d = vel(1)/g*(vel(2) + sqrt(vel(2)^2 + 2*g*pos(2))) + pos(1);
