function [r, dr] = extract_velocity(ref, dt)
% Creates from the position reference 'ref' a position + velocity reference
% as dr(i) = (r(i+1) - r(i-1))/dt^2
%
% Input:
%	ref: the position reference [n, d]
%	dt: the time difference [1]
%
% Output:
%	r: the new position reference [n-2, d]
%	dr: the new velocity reference [n-2, d]

[n, d] = size(ref);
r = zeros(n-2, d);
dr = zeros(n-2, d);

for i = 2:n-1
	r(i-1, :) = ref(i, :);
	dr(i-1, :) = (ref(i+1, :) - ref(i-1, :))./dt/2;
end
