function multiDimDataAnalysis(x, y, subsample)
%
% Plots the columns of x againts the colums of y
%

[n, d] = size(x);
[n, e] = size(y);

x = x(1:subsample:end, :);
y = y(1:subsample:end, :);

figure
for i = 1:e
    for j = 1:d
        subplot(e, d, (i-1)*d+j)
        plot(x(:, j), y(:, i), '+')
    end
end

