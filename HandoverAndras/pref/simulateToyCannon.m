function r = simulateToyCannon(a, s)

x = 0:.01:12;
angleNoise = 1/180*pi;

% Create hills
c = [3 5 6 7];
h = [.3 .3 .4 .2];
scale = [.5 .5 .1 2];
y  = zeros(size(x));
ycontext = 0;

for i = 1:length(c)
    y = y + h(i) * exp(-(x-c(i)).^2/scale(i));
    ycontext = ycontext + h(i) * exp(-(s-c(i)).^2/scale(i));
end
figure, plot(x, y)

% centers for RBF policy
policyc = [2 4 6 8];

for i = 1:length(policyc)
   w(i) = exp(-(s-policyc(i))^2/2);
end

a = a(:);
    
% action and fixed velocity
th = a(1) + w(:)'*a(2:end)/sum(w) + randn(1)*angleNoise;
v = 1;

yproj = x.*tan(th) - 9.81/100 * x.^2/2/(v*cos(th))^2;

ixLand = find(bsxfun(@lt, yproj, y));
ixLand = ixLand(2);
hold on, plot(x(1:ixLand), yproj(1:ixLand), 'r--')
hold on, plot(s, ycontext, 'ko')

r = -sqrt( (x(ixLand) - s).^2 + (yproj(ixLand) - ycontext)^2);

r = r+ 4; % to address 0 prior rewards
