%Condition number


A=[1 0 0; 0 1 0; 0 0 eps]
b=[1 0 0]'

%a
% 2-norm
% norm(A,2) = sqrt(eig_max(A'A))=(since A is diagonal A'A =
% A^2)=sqrt(eig_max(A^2)) = sqrt((eig_max(A))^2) = eig_max(A) = 1.
% For the same reasons, norm(A^-1,2) = 1/epsilon
% => K2 = norm(A,2) * norm(A^-1,2) = 1/epsilon.

% 1-norm / infinity-norm
% Since A is symmetric, the two norms coincide. In particular:
% norm(A,1) = max_j(sum_i(abs(aij))) = 1
% norm(A^-1,1) = max_j(sum_i(abs(aij^-1))) = 1/epsilon
% K1 = Kinf = 1/epsilon

A = eye(3);
epsilon = 1e-6;
A(3,3) = epsilon;


cond(A,2)

cond(A,1)

cond(A,'inf')


%b,c

x_ex= A\b

for (alpha = [1e-6 1e-12])
	deltab_1 = [0 0 alpha]';
	x1 = A \ (b + deltab_1) % x + delta x
	norm(x1 - x_ex, inf) / norm(x_ex, inf) % ||deltax||/||x||
	alpha/epsilon
	
% 	deltab_2 = [alpha 0 0]';
% 	x2 = A \ (b + deltab_2)
% 	norm(x2 - x_ex, inf) / norm(x_ex, inf)
% 	alpha
end





%Iterative methods
A1 = [1 2 2; 
      1 1 1; 
      2 2 1];      
A2 = [2 1 -2; 1 2 1; 2 1 2];   
A3 = [1 2 -2; 1 1 1; 2 2 1];      
A4 = [1 1 -2; 1 2 1; 2 1 2]; 
x=[1 2 3]'

b=A*x

tol=1.e-6
max_it=100

x0 =[0 0 0]'

%Jacobi

D = diag( diag(A1) )
B= inv(D)*(D-A)
g=inv(D)*b

P=D

[x, iter, incr] = stationary_method(B, g, x0, tol, max_it)



rhoBj = max( abs( eig(B) ) )

%Gauss-Seidel
E = -tril(A1, -1)

B=inv(D-E)*(D-E-A)
g=inv(D-E)*b

P= (D-E)

[x, iter, incr] = stationary_method(B, g, x0, tol, max_it)



rhoBj = max( abs( eig(B) ) )







%6.3
n =10;
A =3* eye(n) -2* diag ( ones (n -1 ,1) ,1) - diag ( ones (n -1 ,1) , -1)
b = A * ones (n ,1);

D = diag( diag(A) )
E = -tril(A, -1)

Bgs=inv(D-E)*(D-E-A)
ggs=inv(D-E)*b

Bj=inv(D)*(D-A)
gj=inv(D)*b

pgs=max(abs(eig(Bj)))
pj=max(abs(eig(Bgs)))

x0 = zeros(10,1);
tol = 1e-12;
maxit = 500;

% Jacobi method
[xj, iterj, incrj] = stationary_method(Bj, gj, x0, tol, maxit)

% Gauss-Seidel method
[xgs, itergs, incrgs] = stationary_method(Bgs, ggs, x0, tol, maxit)