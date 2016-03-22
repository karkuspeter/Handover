
ip = '172.31.1.147';
%ip = '127.0.0.1';
tcpipClient = tcpip(ip, 30010,'NetworkRole','Client');

set(tcpipClient,'InputBufferSize',7688);
set(tcpipClient,'Timeout',60);

fopen(tcpipClient);

[context, res] = tcp_getcontext(tcpipClient)

data = tcp_getdata(tcpipClient)

tcp_runexperiment(tcpipClient, context+1)

fclose(tcpipClient);
delete(tcpipClient);
%rawData = fread(tcpipClient,961,'double');
