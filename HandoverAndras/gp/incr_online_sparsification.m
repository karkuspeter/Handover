function [S, update, ix_update] = incr_online_sparsification(dnew, thres, N, kernel_pars, S,verb)
% 
%	Incremental Online Sparsification of training data.
%	Used for selecting the most important training samples 
%	(maximum of N) for incremental GPR.
%
%	Based on: Incremental Online Sparsification for Model Learning
%			  in Real-time Robot Control (Ngyuen-Tuong D. and Peters J., 2011)
%
%	input:
%		dnew: the new data point [1 x n]
%		thres: the distance threshold to include a new sample [1]
%		N: the maximum number of training samples [1]
%		kernel_pars: a callable kernel function with two [1 x n] 
%		S: Data structure with elements
	%		D: dictionary so far [t x n]
	%		a: a list of vectors used for sparsification [t x t-1]
	%		k: a list of vectors used for sparsification [t x t-1]
	%		delta: the distance measure [t]
%		verb: verbosity [1]
%	
%	output:
%		The updated D, Kinv, a, k, delta, with maximum size of N
%		update: if Dictionary is updated 1 , 0 o.w.
%		ix_update: tells, which element of the Dictionary is changed to the new element
	if isfield(S, 'D')
		members = size(S.D, 1);
	else
		members = 0;
	end

	update = 0;

	W = diag(kernel_pars.^-2);


	if members == 0 % init the matrices

		S = {};
		S.D = dnew;
		S.a = [];
		S.k = [];
		S.delta = [];
		update = 1;
		ix_update = 1;

	elseif members == 1 % init both for 1 and 2 the matrices/vectors

		K = 1; 
		k = exp(-.5*(S.D - dnew)*W*(S.D - dnew)');
		a = k/K;

		delta = 1 - k*a;
		
		S.D = [S.D; dnew];
		S.delta = [delta, delta];
		S.a = zeros(2, N-1);
		S.K = zeros(2, N-1);
		S.a(:, 1) = [a; a];
		S.k(:, 1) = [k; k];
		S.Kinvi = zeros(N-1, N-1, 2);
		S.Ki = zeros(N-1, N-1, 2);
		S.Ki(1, 1, 1) = K;
		S.Ki(1, 1, 2) = K;
		S.Kinvi(1, 1, 1) = 1/K;
		S.Kinvi(1, 1, 2) = 1/K;
		S.K = [1, k; k, 1];
		S.Kinv = inv(S.K);

		update = 1;
		ix_update = 2;

	else

		if verb
			disp(' ')
			disp(' ')
		end



		k = exp(-.5 *diag((S.D - repmat(dnew, members, 1))*W*(S.D-repmat(dnew, members, 1))'))';		
		a = S.Kinv*k';
		delta = 1 - k*a(:);

		if verb
			disp(['Distance of new training point ' , num2str(dnew), ' from space spanned by dictionary entries: ',num2str(delta)])
		end

		if delta > thres
			if verb
				disp(['The distance is bigger then threshold (',num2str(delta),'>',num2str(thres),')'])
			end
			if members < N % just inlcude and update
				if verb
					disp(['Adding to the dictionary!'])
				end

				
				for i=1:members
					kmp1v = [k(1:(i-1)), k((i+1):end)];  % k(D^i, d_{m+1}) (vector)
					kimp1 = k(i); % k(d{m+1}, d_i)
					alphai = S.Kinvi(1:members-1, 1:members-1, i)*kmp1v(:);
					gammai = 1 - kmp1v(:)'*alphai(:);

					% Eq. 8
					ainew = [gammai*S.a(i, 1:members-1)' + alphai(:)*alphai(:)'*S.k(i, 1:members-1)' - kimp1*alphai; -alphai'*S.k(i, 1:members-1)'+kimp1]/gammai;
					kinew = [S.k(i, 1:members-1), kimp1]; % Eq. 6
					deltai = 1 - kinew*ainew;

					% Eq. 7
					Kinvnew = 1/gammai*[ gammai*S.Kinvi(1:members-1, 1:members-1, i) + alphai*alphai', -alphai; -alphai', 1];

					% Update data
					S.Ki(1:members, 1:members, i) = [S.Ki(1:members-1, 1:members-1, i), kmp1v(:); kmp1v(:)', 1];
					S.Kinvi(1:members, 1:members, i) = Kinvnew;
					S.a(i, 1:members) = ainew(:)';
					S.k(i, 1:members) = kinew(:)';
					S.delta(i) = deltai;
				end

				% add new training point data
				S.Ki(1:members, 1:members, members+1) = exp(-.5*maha(S.D, S.D, W));
				S.Kinvi(1:members, 1:members, members+1) = inv(S.Ki(1:members, 1:members, members+1));
				S.D = [S.D; dnew];
				S.K = exp(-.5*(maha(S.D, S.D, W)));
				S.Kinv = inv(S.K);
				S.a(members+1, 1:members) = a(:)';
				S.k(members+1, 1:members) = k(:)';
				S.delta = [S.delta, delta];
				ix_update = members+1;
				update = 1;

				if verb 
					disp('Dictionary parameters updated!')
				end

			else % if members == N and delta_new > min(delta) --> replace one

				if delta > min(S.delta)
					if verb
						disp(['New data point is not the least informative --> updating! (delta_min = ', num2str(min(S.delta)),')'])
					end

					ix = find(S.delta == min(S.delta));
					if length(ix) > 1
						ix = ix(1);
					end
					ix_update = ix;
					update = 1;

					S.D(ix, :) = dnew;
					k = exp(-.5 *diag((S.D - repmat(dnew, members, 1))*W*(S.D-repmat(dnew, members, 1))'))';		
					S.K = exp(-.5 * maha(S.D, S.D, W));
					S.Kinv = inv(S.K);

					S.k(ix, :) = [k(1:ix-1), k(ix+1:end)];
					S.a(ix, :) = (S.Kinvi(:, :, ix)*S.k(ix, :)')';
					S.delta(ix) = 1 - S.k(ix, :)*S.a(ix, :)';

					for i = 1:N
						if i ~= ix
							kmp1v = [k(1:i-1), k(i+1:end)];
							kimp1 = k(i);

							if i < ix
								j = ix -1;
							else 
								j = ix;
							end

							Koldjrow = S.Ki(j, :, i);
							r = kmp1v - Koldjrow;

							% Eq. 10
							
							Astar = S.Kinvi(:, :, i) - (S.Kinvi(:, :, i) * r(:)*S.Kinvi(j, :, i))/ (1 + r(:)'* S.Kinvi(j, :, i)');
							A = Astar - (Astar(:, j) * r(:)' * Astar)/(1 + r(:)'*Astar(:, j));

							kinew = [S.k(i, 1:j-1), kimp1, S.k(i, (j+1):end)];

							ainew = A*kinew(:);
							deltai = 1 - kinew*ainew(:);

							% Update data
							S.Kinvi(:, :, i) = A;
							S.Ki(:, :, i) = ...
							[S.Ki(1:j-1, 1:j-1, i), kmp1v(1:j-1)', S.Ki(1:j-1, j+1:end, i);
							kmp1v(1:j-1), 1, kmp1v(j+1:end);
							S.Ki(j+1:end, 1:j-1, i), kmp1v(j+1:end)', S.Ki(j+1:end, j+1:end, i)];
							
							S.a(i, :) = ainew(:)';
							S.k(i, :) = kinew(:)';
							S.delta(i) = deltai;
						end
					end
				
				else % members == N but not good enough
					if verb
						disp(['New delta is higher than threshold, but still the least informative --> no update! (delta_min = ', num2str(min(S.delta)),')'])
					end
					update = 0;
					ix_update = -1;
				end
			end
		else
			ix_update = -1;
			update = 0;
			if verb
				disp('New data is below threshold  -- > update = 0')
			end
		end
	end % members loop
