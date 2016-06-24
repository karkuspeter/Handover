addpath('Matlab_network');

ip = '172.31.1.147';
port = 30000;

%% set up tcp client
tcpClient = tcpip(ip, port, 'NetworkRole','Client');
set(tcpClient,'InputBufferSize',7688);
set(tcpClient,'Timeout',60);

fopen(tcpClient);

%[context, res] = tcp_getcontext(tcpClient)
%data = tcp_getdata(tcpClient)


%% learning
% init REPS, load up if needed

while 1
    % get theta sample from REPS
    theta = [100 -1 5];
    
    % get confirmation from user
    disp('Execute with theta = ');
    disp(theta);
    disp(' ? \n');
    string = input('Continue? [Y]/N ','s');
    if strcmp(string,'Y') || strcmp(string,'y') || isempty(string)
        % can continue
    else
        break;
    end
    
    tcp_runexperiment(tcpClient, theta)

    % get score manually from user
    string = '';
    while 1
        string = input('Input score: ','s');
        score = str2num(string);
        break
    end
    
    % feed it to REPS
    score
    
    % save state in case
end

% evaluate

%% clean up
fclose(tcpClient);
delete(tcpClient);