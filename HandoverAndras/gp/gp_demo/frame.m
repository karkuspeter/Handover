clear all, close all

addpath ../
addpath ~/svnprojects/ClassSystem/Helpers/shadedErrorBar/

x = [-2:.21:10]';
y = sin(x)./x + randn(length(x), 1) * .05;

% figure, plot(x, y, '*')

w = 2.5;
sigf = .2;
sign = .05;

randix = randperm(length(x));
N = 5;
fs = 12;
randNums = randn(length(x), N);


mu = zeros(length(x), 1);
cov = sigf^2 * exp(-.5 * maha(x, x, diag(w.^-2)));
figure,
for j = 1:N
    fapprox(:, j) = cov * randNums(:, j);
end
plot(x, fapprox, 'LineWidth', 1), hold on
[h] = shadedErrorBar(x, mu, 2*diag(cov).^.5, '-k', 1);
    set(h.patch, 'FaceAlpha', .2);
    set(h.mainLine, 'LineWidth', 2);
    set(h.mainLine, 'visible', 'off');
    set(h.edge(1), 'Visible', 'off');
    set(h.edge(2), 'Visible', 'off');
    axis([-2, 10, -.5, .5])
    xlabel('x', 'FontSize', fs), ylabel('y', 'FontSize', fs)
set(gca, 'FontSize', fs)

xtrain = (rand(10, 1) - 1/6) * 12;
ytrain = sin(xtrain)./xtrain + randn(length(xtrain), 1) * .05;

n = 2;
for i = 1: floor(length(xtrain)/n)
  

    K = sigf^2 * exp(-.5 * maha(xtrain(1:(n*i-1)), xtrain(1:(n*i-1)), diag(w.^-2)));
    KG = sigf^2 * exp(-.5 * maha(x, xtrain(1:(n*i-1)), diag(w.^-2)));
    G = sigf^2 * exp(-.5 * maha(x, x, diag(w.^-2)));
    
    mu = KG/(K+sign^2*eye(size(K))) * ytrain(1:(n*i-1));
    cov = G - KG/(K+sign^2*eye(size(K)))*KG';
    sigs = diag(cov).^.5;
    
    
    
    figure, 
    [h] = shadedErrorBar(x, mu, 2*sigs, '-k', 1);
    set(h.patch, 'FaceAlpha', .2);
    set(h.mainLine, 'LineWidth', 2);
    set(h.mainLine, 'visible', 'off');
    set(h.edge(1), 'Visible', 'off');
    set(h.edge(2), 'Visible', 'off');
    hold on,
    h = plot(xtrain(1:n*i-1), ytrain(1:n*i-1), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    
    for j = 1:5
        fapprox(:, j) = mu + cov * randNums(:, j);
    end
    
    plot(x, fapprox, 'LineWidth', 1)
    xlabel('x', 'FontSize', fs), ylabel('y', 'FontSize', fs)
set(gca, 'FontSize', fs)
    
    
end
    