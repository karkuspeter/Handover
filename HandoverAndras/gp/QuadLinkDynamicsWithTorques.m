function [xnew, ffwdTorques] = QuadLinkDynamicsWithTorques(x, action, dt)

x = x(:); action = action(:);

masses = [17.5 17.5 26.25 8.75];
lengths = [0.5 0.5 1 1];
inertias = masses.*(lengths.^2 + 0.0001)./12.0;
friction = [0, 0, 0, 0];
g = 9.81;
sim_dt = dt/100;

res = QuadPendulum_C_ForwardModelWithTorques(x, action, dt, ...
            	masses, lengths, inertias, g, friction, sim_dt);
            
xnew = res(1:8);
ffwdTorques = res(9:12);

