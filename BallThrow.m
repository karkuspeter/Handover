addpath('Matlab_network');

ip = '172.31.1.147';
port = 30000;
z_tresh = 100; %100mm over ground
target_coord = [5000, 5000];

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

    % receive trajectory, will wait until timeout
    % should be in the format [x y z alpha beta theta], mm and rad
    % Z is height, positive up
    tr = tcp_getdata(tcpClient);
    
    % compute coords of hitting the ground (first gets under z treshold)
    first_ind = find(tr(:,3) < z_tresh, 1);
    if ~first_ind
        fprintf('Never went under treshold %f. Closest was %f.', z_tresh, min(tr(:,3)));
    end
    hit_coord = tr(first_ind, :);
    diff2target = hit_coord - target_coord;
        
    % plot trajectory
    figure(1);
    stem3(tr(:,1), tr(:,2), tr(:,3));
    hold on;
    stem3(hit_coord(1), hit_coord(2), hit_coord(3), '-r');
    hold off;
    
    fprintf('Difference to target: [%f, %f], %f', diff2target(1), diff2target(2), norm(diff2target))
    
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