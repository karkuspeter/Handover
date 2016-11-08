clear all, close all

try
    matlabpool(8)
end

addpath('~/Dropbox/workonline/progs/matlab/gpml'); 
startup

sinc = @(x) sin(x)./x;

width = 30;

x = -width/2:.11:width/2;
y = sinc(x);

N = 20;
xsamp = width*(rand(N, 1)-.5);
ysamp = sinc(xsamp) + randn(N, 1)*.1;


[ww, ssigf] = meshgrid(.01:.25:6, .1:.1:3.4);

wwShaped = reshape(ww, [], 1);
ssigfShaped = reshape(ssigf, [], 1);

fShaped = zeros(length(wwShaped), 1);
dfShaped = zeros(length(wwShaped), 2);

parfor i = 1:length(wwShaped);

    loghyp = log([wwShaped(i), ssigfShaped(i), .1])';
    
    [f, df] = hyp_optim_kldiv(loghyp, xsamp, ysamp);
    df = df./exp(loghyp);
    df = df(1:2);
    
    fShaped(i) = f;
    dfShaped(i, :) = df(:)';
end


ff = reshape(fShaped, size(ww));

figure, contour(ww, ssigf, ff, 20), hold on
dfShapedScaled = bsxfun(@rdivide, dfShaped, sqrt(sum(dfShaped.^2, 2)));
quiver(wwShaped, ssigfShaped, dfShapedScaled(:, 1), dfShapedScaled(:, 2));
title('Analytic Gradient (scaled)'), xlabel('length scale'), ylabel('sigf')

figure, contour(ww, ssigf, ff, 20), hold on
quiver(wwShaped, ssigfShaped, dfShaped(:, 1), dfShaped(:, 2));
title('Analyti Gradient (non-scaled)'), xlabel('length scale'), ylabel('sigf')

parfor i = 1:length(wwShaped);

    loghyp = log([wwShaped(i), ssigfShaped(i), .1])';
    
    [f, df] = hyp_optim_kldiv_numerical(loghyp, xsamp, ysamp);
    df = df./exp(loghyp);
    df = df(1:2);
    
    fShaped(i) = f;
    dfShaped(i, :) = df(:)';
end


ff = reshape(fShaped, size(ww));

figure, contour(ww, ssigf, ff, 20), hold on
dfShapedScaled = bsxfun(@rdivide, dfShaped, sqrt(sum(dfShaped.^2, 2)));
quiver(wwShaped, ssigfShaped, dfShapedScaled(:, 1), dfShapedScaled(:, 2));
title('Numerical Gradient (scaled)'), xlabel('length scale'), ylabel('sigf')

figure, contour(ww, ssigf, ff, 20), hold on
quiver(wwShaped, ssigfShaped, dfShaped(:, 1), dfShaped(:, 2));
title('Numerical Gradient (non-scaled)'), xlabel('length scale'), ylabel('sigf')