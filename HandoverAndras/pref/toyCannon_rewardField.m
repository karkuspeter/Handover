clear all, close all

%  r = simulateToyCannon_SimplePolicy(K, x0);
 
 x = 2:.1:10;
 k = 0:.01:pi/2;
 
 [xx, yy] = meshgrid(x, k);
 
 xc = reshape(xx, [], 1);
 yc = reshape(yy, [], 1);
 
 [n,m] = size(xx);
 
 
 for i = 1:length(xc)
     
     rc(i) = simulateToyCannon_SimplePolicy(yc(i), xc(i));
     
 end
 
 rr = reshape(rc, n, m);
 
 figure, contourf(xx, yy, rr, 20);
 

figure, surf(xx, yy, rr, 'EdgeColor', 'none')
xlabel('s')
ylabel('\omega')
zlabel('reward')
title('Toy Cannon task')