function [logHypOpt, fXs] = getFullGPModelGPU(input, target, fevals, checkGradient, gpuoff, initLogHyp, sameWidthParams)
%
% Standard GP regression model learning using squared exponential kernel,
% uses GPU to calculate results.
%
% hyp = exp(logHypInit);
% w_i = hyp(1:end-2, i);
% sigf = hyp(end-1, i);
% sign = hyp(end, i);
%
% k(a, b) = sigf^2 * exp(-.5 * (a-b)' * diag(w_i.^-2) * (a-b));
%
% Input:
%   input: training input [n, d]
%   target: training target [n, e]
%   fevals: maximum number of function evaluations [1]
%       - if negative: maximum number of function evaluation with gradient
%       - if positive: maximum number of line searches
%   checkGradient: checks if we get the same gradient with GPU single as
%   with GPU double and GPML cpu implementation
%   initLogHyp: start the optimization from here
% Output:
%   logHypOpt: the optimal logh hyper parameters [d+2, e]
%   fXs: the final objective function values [e]

if nargin < 4
    checkGradient = 0;
    gpuoff = 0;
    initLogHyp = [];
    sameWidthParams = 0;
elseif nargin < 5
    gpuoff = 0;
    initLogHyp = [];
    sameWidthParams = 0;
elseif nargin < 6
    initLogHyp = [];
    sameWidthParams = 0;
elseif nargin < 7
    sameWidthParams = 0;
end


if checkGradient
    addpath('../../Helpers/gpml'); startup
end

try
    disp('Checking GPU...')
    gpuDevice();
    disp('... GPU detected!')
catch
    gpuoff = 1;
    disp('... GPU not detected, using CPU!')
end

[n, d] = size(input);
[n, e] = size(target);

logHypOpt = zeros(d+2, e);
logHypInit = zeros(d+2, e);

if isempty(initLogHyp)
    hypInitX = std(input);
else
    hypInitX = exp(initLogHyp);
end
logHypInit(1:d, :) = log(repmat(hypInitX(:), 1, e));
for i = 1:e
    logHypInit(d+1, i) = log(std(target(:, i)));
    logHypInit(d+2, i) = log(exp(logHypInit(d+1, i))/10);
end

if ~gpuoff
    reset(gpuDevice);
    input = gpuArray(input);
    target = gpuArray(target);
end

total_runtime = 0;
total_fevals = 0;

if ~sameWidthParams
    for i = 1:e
        optimFunc{i} = @(hyp) hyper_optim_GPUoptim_doubleOnly(hyp, input, target(:, i), 0, 1000, gpuoff);
    end
    
    for i = 1:e
        tic;
        if checkGradient
            [res, fX, lines] = minimizeGPU_checkGrad(logHypInit(:, i), optimFunc{i}, fevals, input, target(:, i));
        else
            [res, fX, lines] = minimizeGPU(logHypInit(:, i), optimFunc{i}, fevals);
        end
        fXs(i) = fX(end);
        logHypOpt(:, i) = res;
        finish = toc;
        total_runtime = total_runtime + finish;
        total_fevals = total_fevals + lines;
        SNR = exp(res(d+1))/exp(res(d+2));
        disp(['Runtime (dimension #', num2str(i), '): ', num2str(finish), ' sec, fevals: ', num2str(lines), ', fevals/sec: ', num2str(lines/finish), ', SNR: ', num2str(SNR)]);
    end
    
    disp(['Total runtime (N/d/e: ', num2str(n), '/',num2str(d), '/', num2str(e), '): ', num2str(total_runtime), ' sec, fevals: ', num2str(total_fevals), ', fevals/sec: ', num2str(total_fevals/total_runtime)]);
else
    logHypInitWidth = mean(logHypInit(1:end-2, :), 2);
    logHypInitNoise = reshape(logHypInit(end-1:end, :), [], 1);
    optimFunc = @(hyp) hyper_optim_GPUoptim_doubleOnly_sameWidth(hyp, input, target, 0, 1000, gpuoff);
    
    [res, fX, lines] = minimizeGPU([logHypInitWidth; logHypInitNoise], optimFunc, fevals);
    
    fXs = fX(end);
    
    logHypOpt = repmat(res(1:d), 1, e);
    logHypOpt(d+1, :) = res(d+1:2:end)';
    logHypOpt(d+2, :) = res(d+2:2:end)';
   
end






