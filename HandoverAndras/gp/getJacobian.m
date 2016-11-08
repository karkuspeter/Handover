function [J, si] = getJacobian(theta, environment, numLink)
% For the Quadlink robot this code calculates the Jacobian at the #numLink-th
% link given the input angle of the robot configuration. 
% Outputs: [jacobian, pos]

lengths = environment.lengths;
si = [0 ; 0];
for i = 1:numLink
   si = si + [sin(sum(theta(1:i))); cos(sum(theta(1:i)))] * lengths(i);
end
% si = si + [sin(sum(theta(1:numLink))); cos(sum(theta(1:numLink)))] * lengths(4) / 2;

J  = zeros(2,4);
for j = 0:(numLink - 1)
    pj = [0;0];
    for i = 1:j
       pj = pj + [sin(sum(theta(1:i))); cos(sum(theta(1:i)))] * lengths(i);
    end
    pj = -(si - pj);
    J([1 , 2], j + 1) = [-pj(2) pj(1)];
end
