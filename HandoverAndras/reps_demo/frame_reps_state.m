clear all, close all


s = [0:.05:5];
w = [0:.05:5];

rewFunc = @(c, a) 3*exp(-.5 * ([a;c] - [3.5; 1.5])' / [.8, .4; .4, 1] * ([a;c] - [3.5; 1.5])) + ...
    4*exp(-.5 * ([a;c] - [1.5; 3.0])' / [1, .5; .5, .8] * ([a;c] - [1.5; 3.0])) ;

[ss,ww] = meshgrid([0:.05:5], [0:.05:5]);

ss = reshape(ss, [], 1);
ww = reshape(ww, [], 1);

for i = 1:length(ss);
    rr(i) = rewFunc(ss(i), ww(i));
end

n = length(s);
ss = reshape(ss, [n, n]);
ww = reshape(ww, [n, n]);
rr = reshape(rr, [n, n]);

fs = 18;

s = s(:);
for i = 1:4
    switch i
        case 1
            contextFunc = @(N) 4 + randn(N, 1) .* .4;
            contextPars = [4 .4];
        case 2
            contextFunc = @(N) 1.5 + randn(N, 1) .* .4;
            contextPars = [1.2 .4];
        case 3
            contextFunc = @(N) 2.5 + randn(N, 1) .* .5;
            contextPars = [2.5 .5];
        case 4
            contextFunc = @(N) 5*rand(N, 1);
    end
    [a, A, cov, rew, eta, theta] = reps_state(2.5, 1, contextFunc, rewFunc, .5, 30, 10, [1; randn(3, 1)*.1]);
    
    figure, subplot(6, 6, [3:6, 9:12, 15:18, 21:24])
    
%     if i < 4
%         p = exp(-.5 * (s - contextPars(1)).^2 / contextPars(2)^2);
%         [h3, h4] = contourf(ss, ww, bsxfun(@times, -rr, p'/2));
%         colormap(gray)
%         set(h4, 'LineColor', 'none')
%     end
resolution = 100;
%     grayBackgroundImage(ss, ww, rr, resolution, .7)
%      if i == 4
%         keyboard
%     end
    
    if i < 4
        p = exp(-.5 * (s - contextPars(1)).^2 / contextPars(2)^2);
        pr = bsxfun(@times,  p', rr);
        grayBackgroundImage([ 0 5],[ 5 0], pr, resolution, max(max(pr))/max(max(rr)))
    else
        p = ones(length(p), 1)*.75;
        pr = bsxfun(@times,  p', rr);
        grayBackgroundImage([ 0 5],[ 5 0], pr, resolution, max(max(pr))/max(max(rr)))
    end

%     if i == 4
%         keyboard
%     end
    
    hold on,
    [h1, h2] = contour(ss, ww, flipud((rr - min(min(rr)))/(max(max(rr)) - min(min(rr)))*resolution), 10);
set(gca, 'XTick', [])
set(gca, 'YTick', [])
    hold on,

    hPolicy = plot(s, 5-a-s*A, 'k', 'LineWidth', 2);
        
    if i < 4
        p = exp(-.5 * (s - contextPars(1)).^2 / contextPars(2)^2);
        subplot(6,6,[27:30, 33:36])
        [hDist] = plot(s, p, 'k');
    else
        p = ones(size(s));
        subplot(6,6,[27:30, 33:36])
        [hDist] = plot(s, p, 'k');
    end
    set(gca, 'YAxisLocation', 'right')
    xlabel('Context $s$', 'interpreter', 'latex', 'FontSize', fs)
    
    
    avgRew = []; 
    M = 50;
    randCs = contextFunc(M);
    
    subplot(6,6,[1:2, 7:8, 13:14, 19:20])
    for j = 1:length(w)
       for k = 1:M
           avgRew(k, j) = rewFunc(randCs(k), w(j));
       end
    end
    
    
    hRDist = plot(w, mean(avgRew,1));  hold on,
%     hRUni = plot(w, mean(rr, 2), 'r');
    ylabel('Action $\omega$', 'interpreter', 'latex', 'FontSize', fs)
    
%     set(ylab, 'Position', [ylabPos(1), ylabPos(2) + 3, ylabPos(3)])
    camroll(-270)
    ylab = get(gca, 'ylabel');
    ylabPos = get(ylab, 'Position');
    set(ylab, 'Position', [2.5, ylabPos(2)+.4, ylabPos(3)])
    set(gca, 'XAxisLocation', 'top')
    
    
%     hLegend = legend([hDist, hRDist, hRUni,hPolicy], {'$\tilde{\mu}(s)$','$p^{\mu}(\mathcal{R}_{sa})$', '$p^{\mathcal{U}}(\mathcal{R}_{sa})$', '$E_{\pi}[\omega | s]$'}, 'interpreter', 'latex', 'Location', 'SouthWest');
    hLegend = legend([hDist, hRDist, hPolicy], {'$\tilde{\mu}(s)$','$p^{\mu}(\mathcal{R}_{sa})$',  '$E_{\pi}[\omega | s]$'}, 'interpreter', 'latex', 'Location', 'SouthWest');
%     set(hLegend, 'box', 'off')
    set(hLegend, 'FontSize', fs);
    set(hLegend, 'FontName', 'Fourier');
        


end

rewFunc2 = @(a) rewFunc(contextFunc(1), a);

[mnopt, covopt, rew, Dkl] = reps(2.5, 1, rewFunc2, .5, 30, 20, 1);

% figure,
%  subplot(6, 6, [1:4, 7:10, 13:16, 19:22])
% [h1, h2] = contourc(ss, ww, rr, 20);
% set(h2, 'LineColor', 'none')
subplot(6, 6, [3:6, 9:12, 15:18, 21:24])
hold on,
hPolicy2 = plot(s, ones(size(s))*(5-mnopt), 'k--', 'LineWidth', 2);

hLegend = legend([hDist, hRDist, hPolicy, hPolicy2], {'$\tilde{\mu}(s)$','$p^{\mu}(\mathcal{R}_{sa})$', '$E_{\pi}[\omega | s]$',  '$E_{\pi}[\omega]$'}, 'interpreter', 'latex', 'Location', 'SouthWest');
%     set(hLegend, 'box', 'off')
    set(hLegend, 'FontSize', fs);
    set(hLegend, 'FontName', 'Fourier');

