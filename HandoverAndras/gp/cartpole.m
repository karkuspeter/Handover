function dy = cartpole(t, y, u)

m1 = .5; m2 = .5; l = .6; b = .1; g = 9.81;
dy = zeros(4, 1);

dy(1) = y(2);
dy(2) = (2*m2*l*y(3)^2*sin(y(4)) + 3*m2*g*sin(y(4))*cos(y(4)) + 4*u - 4*b*y(2))/(4*(m1+m2)-3*m2*cos(y(4))^2);
dy(3) = (-3*m2*l*y(3)^2*sin(y(4))*cos(y(4))-6*(m1+m2)*g*sin(y(4))-6*(u-b*y(2))*cos(y(4)))/(4*l*(m1+m2)-3*m2*l*cos(y(4))^2);
dy(4) = y(3);

