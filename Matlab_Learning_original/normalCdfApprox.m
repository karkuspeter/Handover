function normalCdfApprox

x = [-4:.01:4]';


orders = [ 20 ];
colors = ['b', 'r', 'k', 'g', 'b--', 'r--', 'k--'];
figure,
%% Zelen & Severo
b = [0.2316419,  0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429];
t = 1./(1+b(1) .* x);
result = 1 - normpdf(x).*(b(2)*t + b(3)*t.^2 + b(4)*t.^3 + b(5)*t.^4 + b(6)*t.^5 );

figure, plot(x, result)
%% Taylor
% for j = 1:length(orders)
% 
% sum=x;
% value=x;
%   for i=1:orders(j)
%     
%       value=(value.*x.*x/(2*i+1));
%       sum=sum+value;
%   end
%   result=0.5+(sum./sqrt(2*pi)).*exp(-(x.*x)/2);
%   
%   hold on, plot(x, result, colors(j))
%   
% 
% end