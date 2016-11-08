clear all, close all

addpath('../gpml');
startup
addpath('~/svn_projects/REPS/matlab/ClassSystem/Helpers/')

x = randn(8, 1)*.3; x(7) = x(7)+pi;
xd = zeros(8, 1); xd(7) = pi;
action = zeros(4, 1);
dt = .01;

K = zeros(4, 8);
K(1, 1) = 4000;
K(1, 2) = 100;
K(2, 3) = 4000;
K(2, 4) = 100;
K(3, 5) = 4000;
K(3, 6) = 100;
K(4, 7) = 4000;
K(4, 8) = 100;

K = K/8;

lgp = {{},{},{},{}};
spgp = {};

d = 8;
e = 4;

hypInit = [
    0.4593    0.4760    0.5211    0.4552
    2.2986    1.9908    1.9538    1.5500
    1.1198    0.0015    0.7054   -0.4977
    2.2301    2.9818    2.9580    2.6274
   -0.2019    0.9297   -0.1341    0.6641
    6.8638    6.3205    6.0891    4.3227
   -0.5674   -0.8977   -0.6073   -0.7607
    4.5046    4.5344    4.1563    4.0243
    6.3609    5.7557    5.4583    3.6798
    1.2597    1.1411    1.1369    1.0985];
hypInitOrig = hypInit;

xsave = [];
tsave = [];
tpredsave = [];
usave = [];

for j = 1:10
    x = 2*randn(8, 1).*[.1, .5, .1, .5, .1, .5, .1, .5]'; x(7) = x(7)+pi;
    tic
    for i = 1:500
        e = x-xd;
        e(1:2:8) = mpi2pi(e(1:2:8));
        xext = [sin(x(1:2:8)); cos(x(1:2:8)); x(2:2:8)]';
        ffwdTorquePredict = localGPPredict(lgp, x', 10);
%         ffwdTorquePredict = zeros(1, 4);
%         if ~isempty(spgp)
%             for i = 1:4
%                 k = exp(spgp.loghyp(d+1, i))*exp(-.5*maha(x', spgp.px, diag(exp(spgp.loghyp(1:d, i)))));
%                 ffwdTorquePredict(i) = k*spgp.mf(:, i);
%             end
%         end
        action = -K*e + randn(4, 1)*1 - ffwdTorquePredict';
        tpredsave = [tpredsave; ffwdTorquePredict];
        usave = [usave; action'];
        [xnew, ffwdTorques] = QuadLinkDynamicsWithTorques(x, action, dt);
        xnew(1:2:8) = mpi2pi(xnew(1:2:8));
        if xnew(7) < 0
            xnew(7) = xnew(7)+2*pi;
        end
        ffwdTorques = ffwdTorques + randn(4, 1)*3;
        xsave = [xsave; x'];
        tsave = [tsave; ffwdTorques'];
        xextnew = [sin(xnew(1:2:8)); cos(xnew(1:2:8)); xnew(2:2:8)]';
        lgp = localGPUpdate(lgp, x', ffwdTorques', .1, 50, 0, hypInit);
        x = xnew;
    end
    jani = toc;
    disp(['Simiulation time: ', num2str(jani)])
    disp([num2str(3/jani),' faster than realtime'])
    figure(1), subplot(3,1,1); plot(xsave(:, 1:2:end)); ylabel('q')
    subplot(3,1,2); plot(xsave(:, 2:2:end)); ylabel('dq')
    subplot(3,1,3); plot(usave); ylabel('tau'), drawnow
    figure(2), subplot(2,1,1), plot(tsave); ylabel('ffwdTorques')
    hold on, plot(tpredsave, '--'); 
    subplot(2,1,2), plot(tsave+usave); ylabel('ffwdTorques + fbTorques'), drawnow
    
    % Train sparse gp model
%     if j > 5
% %         keyboard
%         xextsave = [sin(xsave(:, 1:2:8)), cos(xsave(:, 1:2:8)), xsave(:, 2:2:8)];
        subsample = 5;
        m = 300/subsample;
        
        
        if isempty(spgp)
            ix = randperm(size(xsave(1:subsample:end, :), 1));
            initpx = xsave(1:subsample:end, :);
            [spgp.loghyp, spgp.px] = hyper_variationalPseudo_commonInputs(xsave(1:subsample:end, :), tsave(1:subsample:end, :), m, 400, [], carl2ed(hypInit), initpx(ix));
        else
            [spgp.loghyp, spgp.px] = hyper_variationalPseudo_commonInputs(xsave(1:subsample:end, :), tsave(1:subsample:end, :), m, 200, [], carl2ed(hypInit), spgp.px);
        end
        hypInit = ed2carl(spgp.loghyp);
%         spgp.mf = variationalSparsePrecomputations(spgp.loghyp, spgp.px, xsave, tsave, 1);
%     end
    
end
