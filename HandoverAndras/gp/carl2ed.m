function hyped = carl2ed(hypcarl)
% 
% Converts the hyperparameters from 'Rasmussen Style' to 'Snelson Style'
%	log(w, sign, sigf) --> log(b, c, sig)
%

hyped = hypcarl;
hyped(1:end-2, :) = hyped(1:end-2, :)*-2;
hyped(end-1, :) = hyped(end-1, :)*2;
hyped(end, :) = hyped(end, :)*2;
