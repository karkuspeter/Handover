function res = mahaGPU(a, b, Q)
% Mahabolic Distance on gpu
%
% res = (a-b)^T * Q * (a-b), elementwise
%
% input: 
%   a: array of size [n, m];
%   b: array of size [l, m];
%   Q: array of size [m, m];
% ouptut:
%   res: array of size [n, l];

