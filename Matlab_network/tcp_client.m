
%ip = '172.31.1.147';
ip = '127.0.0.1';
tcpipClient = tcpip(ip, 6788,'NetworkRole','Client');

set(tcpipClient,'InputBufferSize',7688);
set(tcpipClient,'Timeout',15);

fopen(tcpipClient);

%msg = read(tcpipClient);
msg = fscanf(tcpipClient);
fprintf('[SERVER]: %s\n', msg);

fprintf(tcpipClient, 'Params\n');
fwrite(3.234);


msg = fscanf(tcpipClient);
fprintf('[SERVER]: %s\n', msg);


fclose(tcpipClient);
delete(tcpipClient);
%rawData = fread(tcpipClient,961,'double');
