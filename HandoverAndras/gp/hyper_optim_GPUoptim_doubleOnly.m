function [f, g] = hyper_optim_GPUoptim_doubleOnly(loghyp, input, target, verb, scaleExtra, gpuoff)

% hyp: the hyperparameters w, sigf, sign
% input: n x d matrix
% target: n x 1 vector
% verb: verbosity

[n, d] = size(input);

w = exp(loghyp(1:d));
sigf = exp(loghyp(d+1));
signs = exp(loghyp(d+2));

% if (sigf/signs) > 2000
%     keyboard
% end

W = diag(w.^-2);
K = sigf*sigf * exp(-.5*maha(input, input, W));
if ~gpuoff
    Ky = K + gpuArray.eye(n)*signs*signs;
else
    Ky = K + eye(n)*signs*signs;
end


L = chol(Ky);

if ~gpuoff
    iKy = (gpuArray.eye(n)/L)/L';
else
    iKy = (eye(n)/L)/L';
end

ldKy = 2*sum(log(diag(L)));

alpha = iKy*target;

f = -(-.5*target'*alpha-.5*ldKy - n/2*log(2*pi)) ;

k = 2; th = 100;
SNR = sigf/signs;
SNRpenalty = -log(factorial(k)) -k*log(th) + (k-1)*log(SNR) - SNR/th;
scale = -length(target)/5 * scaleExtra;

if SNR < 1500
    scale = 0;
end

f = f + scale*SNRpenalty;


if nargout > 1
    dummy = alpha*alpha' - iKy;
    
    if ~gpuoff
        g = gpuArray.zeros(length(loghyp), 1);
    else
        g = zeros(length(loghyp), 1);
    end
    
    
    dSNRpenalty = scale*(-1/th + (k-1)/SNR)*SNR/2;
    
    dsign_snr = -dSNRpenalty;
    dsigf_snr = dSNRpenalty;
    
    % dK / dsigf
    p1 = (2/sigf*K *sigf);
    p1 = single(dummy)*single(p1);
    
    g(end-1, 1) = -.5*trace(double(p1)) + dsigf_snr;
    % dK / dsign
    g(end, 1) =  -.5*trace(2*signs*dummy*signs) + dsign_snr;
    % dK / dwi
    
    for i = 1:d
        if ~gpuoff
            Wloc = gpuArray.zeros(d);
        else
            Wloc = zeros(d);
        end
        
        Wloc(i, i) = 1/w(i)/w(i)/w(i);
        p3 = maha(input, input, Wloc);
        p3 = single((K.* p3)* w(i));
        p3 = single(dummy)*p3;
        g(i, 1) = -.5*double(trace(p3));
    end
    
    if ~gpuoff
        f = (gather(f));
        g = (gather(g));
    end
    if verb
        f, g
    end
    
    
    % keyboard
end
end
