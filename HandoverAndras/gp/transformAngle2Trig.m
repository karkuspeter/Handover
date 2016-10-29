function [mn, sn] = transformAngle2Trig(m, s)
% Transforms angle data to trigonometric functions of the angle. 
% (used for propagating uncertainty in GP trajectory prediction)
% 
% Input:
%	m: mean vector (q, qd, u) [3*d, 1]	
%	s: covariance matrix [3*d, 3*d]
%
% Output:
%	mn: the transformed vector (sin(q), cos(q), qd, u) [4*d, 1]
%	sn: the transformed covariance matrix [4*d, 4*d]

d3 = length(m);
d = d3/3;

[mt, st] = trigaug(m, s, 1:d);

sincoscov = st(3*d+1:end, 3*d+1:end);
sincovdqucov = st(3*d+1:end, d+1:3*d);

sn = zeros(4*d, 4*d);
sn(1:2*d, 1:2*d) = sincoscov;
sn(2*d+1:end, 2*d+1:end) = s(d+1:end, d+1:end)
sn(2*d+1:end, 1:2*d) = sincovdqucov';
sn(1:2*d, 2*d+1:end) = sincovdqucov;

mn =...


