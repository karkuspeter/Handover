function r = simulateToyCannon_SimplePolicy(a, s)
x = 0:.01:12;
angleNoise = 1/180*pi;

c = [3 5 6 7];
h = [.3 .3 .4 .2];
scale = [.5 .5 .1 2];
y  = zeros(size(x));
ycontext = 0;

for i = 1:length(c)
    y = y + h(i) * exp(-(x-c(i)).^2/scale(i));
    ycontext = ycontext + h(i) * exp(-(s-c(i)).^2/scale(i));
end

angleNoise = 1/180*pi;

th = a + randn(1)*angleNoise;

th = max(min(th, 89/180*pi), 1/180*pi);

v = 1;

yproj = x.*tan(th) - 9.81/100 * x.^2/2/(v*cos(th))^2;

ixLand = find(bsxfun(@lt, yproj, y));
ixLand = ixLand(2);

figure, plot(x, y)
hold on, plot(x(1:ixLand), yproj(1:ixLand), 'r--')
hold on, plot(s, ycontext, 'ko')
xlabel('Context (distance) [m]')
ylabel('Height [m]')
title('Toy Cannon problem')
legend('Mountain surface', 'Projectile Trajectory', 'Target')

r = -sqrt( (x(ixLand) - s).^2 + (yproj(ixLand) - ycontext)^2);
r = r+ 4;
