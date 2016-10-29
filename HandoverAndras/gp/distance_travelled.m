function d = distance_travelled(th, v, y0);

g = 9.81;
d = v*cos(th)/g*(v*sin(th) + sqrt((v*sin(th))^2 + 2*g*y0));
