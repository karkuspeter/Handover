function [f, g] = hyper_optim_GPUoptim_SqExp_Sq(loghyp, input, target, verb, maximize, gpuoff, numerical)

% hyp: the hyperparameters w, sigf, sign
% input: n x d matrix
% target: n x 1 vector
% verb: verbosity
% maximize: 1 if you use maximizatoin function, -1 if you use minimiaztion
% gpuoff: use only cpu

if nargin < 4
    verb = 0;
    maximize = -1;
    gpuoff = 0;
    numerical = 0;
elseif nargin < 5
    maximize = -1;
    gpuoff = 0;
    numerical = 0;
elseif nargin < 6
    gpuoff = 0;
    numerical = 0;
elseif nargin < 7
    numerical = 0;
end

[n, d] = size(input);

w = exp(loghyp(1:d));
a = exp(loghyp(d+1));
b = exp(loghyp(d+2));
signs = exp(loghyp(d+3));

W = diag(w.^-2);
Ka = a * exp(-.5*maha(input, input, W));
Kb = b * input*input';
if ~gpuoff
    Ky = Ka + Kb + parallel.gpu.GPUArray.eye(n)*signs*signs;
else
    Ky = Ka + Kb + eye(n)*signs*signs;
end

L = chol(Ky);

if ~gpuoff
    iKy = (parallel.gpu.GPUArray.eye(n)/L)/L';
else
    iKy = (eye(n)/L)/L';
end

ldKy = 2*sum(log(diag(L)));

alpha = iKy*target;

f = -(-.5*target'*alpha - .5*ldKy - n/2*log(2*pi)) ;

% k = 2; th = 100;
% SNR = sigf/signs;
% SNRpenalty = -log(factorial(k)) -k*log(th) + (k-1)*log(SNR) - SNR/th;
% scale = -length(target)/5 * scaleExtra;
% 
% if SNR < 1500
%     scale = 0;
% end
% 
% f = f + scale*SNRpenalty;


if and(nargout > 1, ~numerical)
    dummy = alpha*alpha' - iKy;
    
    if ~gpuoff
        g = parallel.gpu.GPUArray.zeros(length(loghyp), 1);
    else
        g = zeros(length(loghyp), 1);
    end
    
    
%     dSNRpenalty = scale*(-1/th + (k-1)/SNR)*SNR/2;
%     
%     dsign_snr = -dSNRpenalty;
%     dsigf_snr = dSNRpenalty;
    
    % dK / da
    g(end-2, 1) = .5*trace(dummy * Ka/a );
    % dK / db
    g(end-1, 1) = .5*trace(dummy * Kb/b );
    % dK / dsign
    g(end, 1) =  .5*trace(dummy*2*signs);
    % dK / dwi
    
    for i = 1:d
        if ~gpuoff
            Wloc = parallel.gpu.GPUArray.zeros(d);
        else
            Wloc = zeros(d);
        end
        Wloc(i, i) = 1/w(i)/w(i)/w(i);
        
        g(i, 1) = .5*trace(dummy* (Ka.*maha(input, input, Wloc)));

    end
    
    if ~gpuoff
        f = (gather(f));
        g = (gather(g));
    end
    
    g = -g.*exp(loghyp(:));
    
   
    if ~maximize
        f = -f;
        g = -g;
    end
    if verb
        f, g
    end
    
elseif and(nargout > 1, numerical)
    
    objfun = @(lhyp) hyper_optim_GPUoptim_SqExp_Sq(lhyp, input, target, verb, maximize, 1);
    g = numericalGradient(objfun, loghyp(:));
    
  
    if ~maximize
        f = -f;
        g = -g;
    end
    
     if verb
        f, g
    end
end
end
