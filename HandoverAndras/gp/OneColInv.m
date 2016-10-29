function Bnew = OneColInv(B,X0,col,posX0,add)
% Written by Emtiyaz, CS, UBC
% Updated on Feb 21, 2008
% X0 is a matrix of size nxp,
% B is inverse of X0 with columns indexed by col i.e. X0(:,col)
% we want to add/remove a new column posX0 and update the inverse. 
% This program implements an update using Matrix-inversion-lemma.
% First add colum to the end, Xnew = [X v]
% let B = inv(X'X) available to us, then using Matrix-inversion-lemma
% and schur complement, we get
% Bnew = inv(Xnew'Xnew) = [F11inv, -B*X'*v*F22inv; -F22inv*B*X'*v F22inv]
% where F11inv = B + B*X'*v*F22inv*v'*X*B
% and F22inv = v'*v - u'*X*B*X'*v
% at the end permute the column and rows, as XnewActual = Xnew*P
% where P is a permutation matrix, and hence
% inv(XnewActual'*XnewActual) = P*(Xnew'*Xnew)*P'

%find the submatrix with col
X = X0(:,col);

if add == 1
  % find the position in X created with col
  pos = find(sort([col; posX0])==posX0);
  % find the column to be inserted
  v = X0(:,posX0);
  % compute the inverse as if adding a colum to the end
  u1 = X'*v;
  u2 = B*u1;
  F22inv = 1/(v'*v - u1'*u2);
  u3 = F22inv*u2;
  F11inv = B + F22inv*u2*u2';
  Bnew = [F11inv -u3; -u3' F22inv];
  % permute to get the matrix corresponding to original X
  Bnew = Bnew(:,[1:pos-1 end pos:end-1]);
  Bnew = Bnew([1:pos-1 end pos:end-1],:);
else
  % find the position in new matrix X 
  pos = find(col==posX0);
  % permute to bring the column at the end in X
  B = B([1:pos-1 pos+1:end pos],:);
  B = B(:,[1:pos-1 pos+1:end pos]);
  %update the inverse by removing the last column
  F11inv =B(1:end-1,1:end-1);
  F22inv = B(end,end);
  u3 = -B(1:end-1,end);
  u2 = u3/F22inv;
  Bnew = F11inv - u2*u2'*F22inv;
end


