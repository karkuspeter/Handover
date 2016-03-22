 tcpipServer = tcpip('0.0.0.0',55000,'NetworkRole','Server');
 
set(tcpipServer,'OutputBufferSize',7688);
fopen(tcpipServer);

%write(tcpipServer, 'Gyurcsany a hibas\n');
fprintf(tcpipServer, 'Gyurcsany a hibas');

msg = fscanf(tcpipServer)


fprintf(tcpipServer, 'csak %s ?', msg);


fclose(tcpipServer);
delete(tcpipServer);
%rawData = fread(tcpipClient,961,'double');
