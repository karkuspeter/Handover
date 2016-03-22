function [ data, res ] = tcp_getdata( tcpipClient )
%Receive matrix. Returns 0 if unsuccessful

%fprintf(tcpipClient, 'Context\n');
msg = fscanf(tcpipClient);
if (~strcmp(msg(1:end-1), 'Data'))
    disp(msg);
    return
end
disp('getting data size');


msg = fscanf(tcpipClient);
rows = str2num(msg);
msg = fscanf(tcpipClient);
cols = str2num(msg);

rows, cols
data = zeros(rows, cols);

for i=1:rows
    for j=1:cols
        msg = fscanf(tcpipClient);
        data(i, j) = str2double(msg);
    end
end

msg = fscanf(tcpipClient);
if (~strcmp(msg(1:end-1), 'End'))
    disp(''+msg);
    data = 0;
    return
end

res =1;
end

