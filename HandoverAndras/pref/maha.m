function K = maha ( a , b , Q )
%
% Squared Mahalanobis d i s t a n c e ( a−b ) *Q * ( a−b ) ' ; v e c t o r s a r e row−v e c t o r s
% a , b m a t r i c e s c o n t a i n i n g n l e n g t h d row v e c t o r s , d by n
%Qweight matrix , d by d , d e f a u l t eye ( d )
% Ksquared d i s t a n c e s , n by n
if nargin == 2
% assume i d e n t i t y Q
K = bsxfun ( @plus , sum ( a .* a , 2 ) ,sum ( b .* b , 2 )' ) -2* a * b' ;
else
aQ = a * Q ; K = bsxfun ( @plus , sum ( aQ .* a , 2 ) ,sum ( b * Q .* b , 2 )' ) -2* aQ * b' ;

end

