function [ context, res ] = tcp_getcontext( tcpipClient )
%Receive context parameters. Returns 0 if unsuccessful
param_count = 3;

context = zeros(param_count,1);
res = 0;
fprintf(tcpipClient, 'Context\n');
msg = fscanf(tcpipClient);
if (~strcmp(msg(1:end-1), 'Context'))
    disp(msg);
    return
end
disp('getting context');

for i=1:param_count
    msg = fscanf(tcpipClient);
    context(i) = str2double(msg);
end

msg = fscanf(tcpipClient);
if (~strcmp(msg(1:end-1), 'End'))
    disp(''+msg);
    context = 0;
    return
end

res =1;
end

