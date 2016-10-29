function [mnopt, covopt, rew, Dkl] = reps(mn, sigs, rewFunc, epsilon, samples, episodes, eta)

options = optimset('Algorithm','active-set');
options = optimset(options, 'GradObj','on');
options = optimset(options, 'Display', 'off');
S = [];
r = [];
N = samples;

dim = length(mn);
cov = diag(sigs.^2);

for e = 1:episodes

    Snew = mvnrnd(mn(:), cov, N);
    
    
	%q = exp(-.5*diag((Snew - repmat(mn', N, 1))*inv(cn)*(Snew - repmat(mn', N, 1))'));
	%q = q/sum(q);
	q = ones(N, 1)/N;
    
    rnew = zeros(N, 1);
    for i = 1:N
        rnew(i) = rewFunc(Snew(i, :));
    end
    
	%S = [S; Snew];
	%r = [r; rnew];
    
	
	S = Snew;
	r = rnew;
    
    ixBad = find(isnan(r));
    ixGood = setdiff(1:length(r), ixBad);
    
    if ~isempty(ixBad)
        warning(['There are ', num2str(length(ixBad)), ' NaN values in the reward! These samples will be discarded'])
    end
    
    S = S(ixGood, :);
    r = r(ixGood);
    

	% With Gradient
	objfun = @(eta) reps_dual(eta, r, q, epsilon);
%     keyboard
    try
        eta_prev = eta;
        eta = fmincon(objfun, eta, -1, -.01, [], [], [], [], [], options);
        if isnan(eta)
            disp('rewards')
            r
            error('NaN eta')
        end
    catch err
        if mean(r)/std(r) > 100
            disp('Optimal solution found!')
            break
        end
    end
        
   
    

	p = exp(r/eta)/sum(exp(r/eta));
	Dkl(e) = sum(p.*log(p./q(ixGood)));
% 	iDkl = sum(q.*log(q./p));
	mn = S'*p;
	cov = bsxfun(@minus, S, mn')'*bsxfun(@times, bsxfun(@minus, S, mn'), p(:));
    
    if any(isnan(mn))
        keyboard
    end

    rew(e) = mean(rnew);
% 	rexp = r(:)'*p;
% 	rvar = (r(:)-rexp)'*((r(:)-rexp).*p);
% 
% 	savem = [savem; mn', diag(cn)'.^.5, cn(2, 1)/prod(diag(cn).^.5), eta, rexp, Dkl, iDkl, rvar];
    disp(['Episode #', num2str(e), ', mean reward: ', num2str(rewFunc(mn(:)'))]);
end
mnopt = mn;
covopt = cov;