function [ M , S , V ] = gpPred ( logh , input , target , m , s )
%
% Compute j o i n t GP p r e d i c t i o n s ( m u l t i v a r i a t e t a r g e t s ) f o r an
% u n c e r t a i n t e s t in pu t x s t a r  ̃ N(m, S )
%
%Ddimension o f t r a i n i n g i n p u t s
% Edimension o f t r a i n i n g t a r g e t
% ns i z e of t r a i n i n g s e t
%
% input arguments :
%
% logh: E * (D+2) by 1 v e c t o r log−hyper−parameters :
%[D log−length−s c a l e s , log−s i g n a l−s t d . dev , log−noise−s t d . dev ]
% input: n by D matrix o f t r a i n i n g i n p u t s
% target: n by E matrix o f t r a i n i n g t a r g e t s
% m: D by 1 mean v e c t o r o f t h e t e s t in pu t
% s: D by D c o v a r i a n c e matrix o f t h e t e s t i np ut
%
% returns :
%
%M: E by 1 v e c t o r , mean o f p r e d i c t i v e d i s t r i b u t i o n  E {h , x }[ h ( x ) |m, s ]eq . ( 2 . 4 3 )
% S: E by E matrix , c o v a r i a n c e p r e d i c t i v e d i s t r i b u t i o n %cov {h , x }[ h ( x ) |m, s ]eq . ( 2 . 4 4 )
% V: D by E c o v a r i a n c e between i n p u t s and outputs
%cov {h , x }[ x , h ( x ) |m, s ]eq . ( 2 . 7 0 )
%
% incorporates :
% a ) model u n c e r t a i n t y about t h e underlying f u n c t i o n ( i n p r e d i c t i o n )
%
%
% Copyright (C) 2008−2010 by Marc P e t e r D e i s e n r o t h and C a r l Edward Rasmussen
% l a s t m o d i f i c a t i o n : 2010−08−30
persistent K iK old_logh ;
% cached v a r i a b l e s
[ n , D ] = size ( input ) ;
% s i z e o f t r a i n i n g s e t and dimension o f in pu t space
[ n , E ] = size ( target ) ;
% s i z e o f t r a i n i n g s e t and t a r g e t dimension
logh = reshape ( logh , D +2 , E )' ;
% un−v e c t o r i z e log−hyper parameters
% i f re−computation n e c e s s a r y : compute K, inv (K) ; o t h e r w i s e use cached ones
if numel ( logh ) ~= numel ( old_logh ) || isempty ( iK ) ...
|| sum ( any ( logh ~= old_logh ) ) || numel ( iK ) ~= E * n^2
old_logh = logh ;
iK = zeros ( n , n , E ) ; K = zeros ( n , n , E ) ;
for i =1: E
inp = bsxfun ( @rdivide , input , exp ( logh ( i , 1 : D ) ) ) ;
K ( : , : , i ) = exp ( 2 * logh ( i , D +1)-maha ( inp , inp ) /2) ;
% k e r n e l matrix K ; n−by−n
L = chol ( K ( : , : , i ) +exp ( 2 * logh ( i , D +2) ) * eye ( n ) )' ;
iK ( : , : , i ) = L' \ ( L\eye ( n ) ) ;
% i n v e r s e k e r n e l matrix inv (K) ; n−by−n
end
end
% memory a l l o c a t i o n
M = zeros ( E , 1 ) ; V = zeros ( D , E ) ; S = zeros ( E ) ;
log_k = zeros ( n , E ) ; beta = zeros ( n , E ) ;

%
%
% steps
% 1 ) compute p r e d i c t e d mean and c o v a r i a n c e between i np ut and p r e d i c t i o n
% 2 ) p r e d i c t i v e c o v a r i a n c e matrix
% 2a ) non−c e n t r a l moments
% 2b ) c e n t r a l moments
inp = bsxfun ( @minus , input , m' ) ; % s u b t r a c t mean o f t e s t i np ut from t r a i n i n g i np ut
% 1 ) compute p r e d i c t e d mean and c o v a r i a n c e between i np ut and p r e d i c t i o n
for i =1:E
% f o r a l l t a r g e t dimensions
beta ( : , i ) = ( K ( : , : , i ) +exp ( 2 * logh ( i , D +2) ) * eye ( n ) ) \target ( : , i ) ;
% K\y ; n−by−1
iLambda = diag ( exp(-2 * logh ( i , 1 : D ) ) ) ;
% i n v e r s e squared length−s c a l e s ; D−by−D
R = s+diag ( exp ( 2 * logh ( i , 1 : D ) ) ) ;
% D−by−D
iR = iLambda * ( eye ( D ) - ( eye ( D ) +s * iLambda ) \( s * iLambda ) ) ;
% Kailath inverse
T = inp * iR ;
% n−by−D
c = exp ( 2 * logh ( i , D +1) ) / sqrt ( det ( R ) ) * exp ( sum ( logh ( i , 1 : D ) ) ) ;
% scalar
q = c * exp(-sum ( T .* inp , 2 ) /2) ;
% eq . ( 2 . 3 6 ) ; n−by−1
qb = q .* beta ( : , i ) ;
% n−by−1
M ( i ) = sum ( qb ) ;
% p r e d i c t e d mean , eq . ( 2 . 3 4 ) ; s c a l a r
V ( : , i ) = s * T' * qb ;
% input−output cov . , eq . ( 2 . 7 0 ) ; D−by−1
v = bsxfun ( @rdivide , inp , exp ( logh ( i , 1 : D ) ) ) ;
% ( X− * s q r t ( iLambda ) ; n−by−m)D
log_k ( : , i ) = 2 * logh ( i , D +1)-sum ( v .* v , 2 ) / 2 ;
% precomputation f o r 2 ) ; n−by−1
end
% 2 ) p r e d i c t i v e c o v a r i a n c e matrix ( symmetric )
% 2a ) non−c e n t r a l moments
for i =1: E
% f o r a l l t a r g e t dimensions
Zeta_i = bsxfun ( @rdivide , inp , exp ( 2 * logh ( i , 1 : D ) ) ) ;
% n−by− D
for j =1: i
Zeta_j = bsxfun ( @rdivide , inp , exp ( 2 * logh ( j , 1 : D ) ) ) ;
R = s * diag ( exp(-2 * logh ( i , 1 : D ) ) +exp(-2 * logh ( j , 1 : D ) ) ) +eye ( D ) ;
t = 1 ./ sqrt ( det ( R ) ) ;
% n−by−D% D−by−D
% scalar
% e f f i c i e n t i m p l e n t a t i o n o f eq . ( 2 . 5 3 ) ; n−by−n
Q = t * exp ( bsxfun ( @plus , log_k ( : , i ) , log_k ( : , j )' ) +maha ( Zeta_i ,-Zeta_j , R\s /2) ) ;
A = beta ( : , i ) * beta ( : , j )' ;
% n−by−n
if i==j % i n c o r p o r a t e model u n c e r t a i n t y ( d i a g o n a l o f S only )
A = A - iK ( : , : , i ) ; % req . f o r E x [ v ar h [ h ( x ) | x ] |m, s ] , eq . ( 2 . 4 1 )
  end 
A = A.*Q; % n−by−n
S ( i , j ) = sum ( sum ( A ) ) ; % i f i == j : 1 s t term i n eq . ( 2 . 4 1 ) e l s e eq . ( 2 . 5 0 )
S(j , i) = S(i , j) ; % copy e n t r i e s
end
% add s i g n a l v a r i a n c e t o diagonal , completes model u n c e r t a i n t y , eq . ( 2 . 4 1 )
S ( i , i ) = S ( i , i ) + exp ( 2 * logh ( i , D +1) ) ;
end
% 2b ) c e n t r a l i z e moments
S = S - M*M' ;
% c e n t r a l i z e moments , eq . ( 2 . 4 1 ) , ( 2 . 4 5 ) ; E−by−E
%− − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
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
%− − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − − −−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
% References
%
% Marc P e t e r D e i s e n r o t h :
% E f f i c i e n t Reinforcement Learning using Gaussian P r o c e s s e s
% PhD Thesis , K a r l s r u h e I n s t i t u t e o f Technology

