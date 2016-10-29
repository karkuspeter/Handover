function [f, g] = hyper_optim_CPUoptim2(loghyp, input, target, verb)

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
Ky = K + eye(n)*signs*signs;

if 0
    go = 1;
    sig_sq_curr = signs*signs;
    while go
        try
            L = chol(Ky);
            go = 0;
        catch
            disp(['chol(Ky) crashing, precision: ', class(gather(Ky)),', SNR: ', num2str(sigf/signs), ', sigf/sign: ', num2str([sigf, signs])])
            sig_sq_curr = sig_sq_curr * 10;
            Ky = Ky + eye(n)*sig_sq_curr;
        end
    end
else
    L = chol(Ky);
end



% if 0
%     iKy = parallel.gpu.GPUArray.eye(n)/L/L';
% else
iKy = single(Ky);
iKy = inv(iKy);
% end
% target = single(target);

ldKy = 2*sum(log(diag(L)));

alpha = iKy*target;

f = -(-.5*target'*alpha-.5*ldKy - n/2*log(2*pi)) ;

k = 2; th = 100;
SNR = sigf/signs;
SNRpenalty = -log(factorial(k)) -k*log(th) + (k-1)*log(SNR) - SNR/th;
scale = -length(target)/5;

if SNR < 500
    scale = 0;
end

f = f + scale*SNRpenalty;


if nargout > 1
    dummy = single(alpha*alpha' - iKy);
    
    g = zeros(length(loghyp), 1);
    
    dSNRpenalty = scale*(-1/th + (k-1)/SNR)*SNR/2;
    
    dsign_snr = -dSNRpenalty;
    dsigf_snr = dSNRpenalty;
    
    % dK / dsigf
    p1 = dummy;
    p1 = (2/sigf*K *sigf);
    p1 = dummy*single(p1);
    
    g(end-1, 1) = -.5*trace(p1) + dsigf_snr;
    % dK / dsign
    g(end, 1) =  -.5*trace(2*signs*dummy*signs) + dsign_snr;
    % dK / dwi
    
    for i = 1:d
        Wloc = zeros(d); Wloc(i, i) = 1/w(i)/w(i)/w(i);

        p3 = maha(input, input, Wloc);
        p3 = single((K.* p3)* w(i));
        p3 = dummy*p3;
        g(i, 1) = -.5*trace(p3);
    end
    
    f = double(f);
    g = double(g);
    if verb
        f, g
    end
    
    
    % keyboard
end
end
