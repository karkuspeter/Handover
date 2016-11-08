function HandoverAndras(send_polpar)

addpath('../Matlab_network');

ip = '172.31.1.147';
port = 30000;

tcpClient = tcpip(ip, port, 'NetworkRole','Client');
set(tcpClient,'InputBufferSize',7688);
set(tcpClient,'Timeout',30);

fopen(tcpClient);

tcp_runexperiment(tcpClient, send_polpar);

fclose(tcpClient);
delete(tcpClient);
