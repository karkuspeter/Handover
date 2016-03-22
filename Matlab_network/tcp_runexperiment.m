function [ res ] = tcp_runexperiment( tcpipClient, params )
%Send run experiment command to KUKA with params

fprintf(tcpipClient, 'Run');

fprintf( '%d\n', length(params));

fprintf(tcpipClient, '%d\n', length(params));
for i=1:length(params)
    fprintf(tcpipClient, '%3.8f\n', params(i));
    fprintf('%3.8f\n', params(i));
end
fprintf(tcpipClient, 'End\n');
res = 1;

% msg = fscanf(tcpipClient);
% if (strcmp(msg, 'ok\n'))
%     disp('reply received');
%     res = 1;
% else
%     disp('no reply');
%     res = 0;
% end
end

