addpath('Matlab_network');

ip = '172.31.1.147';
port = 30000;
%z_tresh = 100; %100mm over ground
target_coord = [5000, 5000];
rigid_body_start_dex = 3;

%% set up tcp client
tcpClient = tcpip(ip, port, 'NetworkRole','Client');
set(tcpClient,'InputBufferSize',7688);
set(tcpClient,'Timeout',30);

fopen(tcpClient);

%[context, res] = tcp_getcontext(tcpClient)
%data = tcp_getdata(tcpClient)


%% learning
% init REPS, load up if needed

while 1
    % get theta sample from REPS
    theta = [-20 0.9];
    
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
%     first_ind = find(tr(:,3) < z_tresh, 1);
%     if ~first_ind
%         fprintf('Never went under treshold %f. Closest was %f.', z_tresh, min(tr(:,3)));
%     end
%     hit_coord = tr(first_ind, :);
%     diff2target = hit_coord - target_coord;

    recsize = size(tr)
    if isempty(tr) || size(tr,2) < 3
        tr = [9 9 9 9; 10 10 10 10];
    elseif rigid_body_start_dex+3 > size(tr,2)
        rigid_body_start_dex = 1;
    end
    final_coord = tr(end,[rigid_body_start_dex rigid_body_start_dex+1]);
    diff2target = final_coord - target_coord;
    
    % plot trajectory
    %figure(1);
    %stem3(tr(:,rigid_body_start_dex), tr(:,rigid_body_start_dex+1), tr(:,rigid_body_start_dex+2));
    %hold on;
    %stem3(hit_coord(1), hit_coord(2), hit_coord(3), '-r');
    %hold off;
    
    fprintf('Difference to target: [%f, %f], %f', diff2target(1), diff2target(2), norm(diff2target))
    
    % get score manually from user
%     string = '';
%     while 1
%         string = input('Input score: ','s');
%         score = str2num(string);
%         break
%     end
%     
    % feed it to REPS
    score = norm(diff2target);
    
    % save state in case
end

% evaluate

%% clean up
fclose(tcpClient);
delete(tcpClient);