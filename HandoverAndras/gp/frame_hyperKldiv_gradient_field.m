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

N = 100;
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

figure, contour(ww, ssigf, ff, 30), hold on
dfShapedScaled = bsxfun(@rdivide, dfShaped, sqrt(sum(dfShaped.^2, 2)));
quiver(wwShaped, ssigfShaped, dfShapedScaled(:, 1), dfShapedScaled(:, 2));
title('KLdiv Gradient (scaled)'), xlabel('length scale'), ylabel('sigf'), colorbar

figure, contour(ww, ssigf, ff, 30), hold on
quiver(wwShaped, ssigfShaped, dfShaped(:, 1), dfShaped(:, 2));
title('KLdiv Gradient (non-scaled)'), xlabel('length scale'), ylabel('sigf'), colorbar

parfor i = 1:length(wwShaped);

    loghyp = log([wwShaped(i), ssigfShaped(i), .1])';
    
    [f, df] = hyper_optim_GPUoptim_doubleOnly(loghyp, xsamp, ysamp, 0, 0, 1);
    df = df./exp(loghyp);
    df = df(1:2);
    
    fShaped(i) = f;
    dfShaped(i, :) = df(:)';
end


ff = reshape(fShaped, size(ww));

figure, contour(ww, ssigf, ff, 30), hold on
dfShapedScaled = bsxfun(@rdivide, dfShaped, sqrt(sum(dfShaped.^2, 2)));
quiver(wwShaped, ssigfShaped, dfShapedScaled(:, 1), dfShapedScaled(:, 2));
title('Standard GP Gradient (scaled)'), xlabel('length scale'), ylabel('sigf'), colorbar

figure, contour(ww, ssigf, ff, 30), hold on
quiver(wwShaped, ssigfShaped, dfShaped(:, 1), dfShaped(:, 2));
title('Standard GP Gradient (non-scaled)'), xlabel('length scale'), ylabel('sigf'), colorbar