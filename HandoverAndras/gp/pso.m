function parOpt = pso(range, n, f, steps, omega, extraConstr)
%
% Particle Swarm Optimization
%
% Input:
%	range: the range of parameters [d, 2]
%	n: number of particles [1]
%	f: objective function [func]
%	steps: the number of steps [1]
%		if negative the problem is minimization
%		if positive the problem is maximization
%   omega: moment parameters \in (0, 1] [3]
%	extraConstr: constraining function with input: particle [list of func]
if nargin < 6
	extraConstr = [];
end

d = size(range, 1);
swarm = zeros(n, d);
swarmBest = zeros(n, d);
objFun = zeros(n, 1);
objFunBest = zeros(n, 1);
absRange = (range(:, 2) - range(:, 1))';

counter = 0;
for i = 1:n
    go = 1;
    
    while go
        swarm(i, :) = rand(1, d) .* (range(:, 2) - range(:, 1))' + range(:, 1)';
        swarmBest(i, :) = swarm(i, :);
        try
            objFun(i) = f(swarm(i, :)');
            objFunBest(i) = objFun(i);
            go = 0;
        catch
            counter = counter + 1;
            [i counter]
        end
    end
end



if steps > 0
	[bestObj, ixBest] = max(objFun);
	[worstObj, ixWorst] = min(objFun);
else
	[bestObj, ixBest] = min(objFun);
	[worstObj, ixWorst] = max(objFun);
end
bestParticle = swarm(ixBest, :);
worstParticle = swarm(ixWorst, :);

disp(['Iteration #0, best/mean/worst obj: ', num2str(bestObj), '/', ...
		num2str(mean(objFun)), '/', num2str(worstObj)]);

vel = bsxfun(@times, 2*(rand(n, d) - .5), absRange);

for i = 1:abs(steps)
	for j = 1:n
	
		acceptable = 0;
        
		while ~acceptable

            rp = rand(1, d);
			rg = rand(1, d);
			newVel = vel(j, :)*omega(1) + omega(2)*rp.*(swarmBest(j, :) - swarm(j, :)) + ...
					omega(3) * rg.*(bestParticle - swarm(j, :));
			
            maxNewVel = range(:, 2)' - swarm(j, :);
            minNewVel = range(:, 1)' - swarm(j, :);
            newVel = min(max(minNewVel, newVel), maxNewVel);
                

            if ~isempty(extraConstr)
                for k = 1:length(extraConstr)
                    acceptableExtra(k) = extraConstr{k}(swarm(j, :) + newVel);
                end
                acceptable = all(acceptableExtra);
            else
                acceptable = 1;
            end
		end
		vel(j, :) = newVel;
		
        try
            swarm(j, :) = swarm(j, :) + vel(j, :);
    		objFun(j) = f(swarm(j, :));
        catch
        end
		
		if steps > 0
			if objFun(j) > objFunBest(j)
				objFunBest(j) = objFun(j);
				swarmBest(j, :) = swarm(j, :);
				if objFunBest(j) > bestObj
					bestObj = objFunBest(j);
					bestParticle = swarmBest(j, :);
				end
			end
		else
			if objFun(j) < objFunBest(j)
				objFunBest(j) = objFun(j);
				swarmBest(j, :) = swarm(j, :);
				if objFunBest(j) < bestObj
					bestObj = objFunBest(j);
					bestParticle = swarmBest(j, :);
				end
			end
		end
	end
	
	if steps > 0
		[worstObj, ixWorst] = min(objFun);
	else
		[worstObj, ixWorst] = max(objFun);
	end
	worstParticle = swarm(ixWorst, :);
	disp(['Iteration #', num2str(i), ', best/mean/worst obj: ', num2str(bestObj), '/', ...
		num2str(mean(objFun)), '/', num2str(worstObj)]);
end

parOpt = bestParticle;
	


	
	
	