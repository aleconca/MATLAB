%5.1 LU 
A=[1 2 1 1; 1 4 0 2; 2 10 4 0; 1 0 2 2]
b=[3 3 10 1]'

[L,U,P]=lu(A)
y=L\(P*b)
x=U\y

%spy(A*P)

%detPA= det(P)*det(A)= det(L)*det(U)

detA= (prod(diag(L))*prod(diag(U)))/det(P)

det(A)





%5.2 Chol

A=[10 0 3 0; 0 5 0 -2;3 0 5 0; 0 -2 0 2]
b=[2 2 2 2]'

eig(A)

R = chol(A)

y = forward_substitution(R',b)  
x = backward_substitution(R,y)

A\b

detA= prod(diag(R))^2;

%5.3 

a=2*ones(10,1)
b=1*ones(10,1)

A = spdiags([b a b], [-1 0 1], 10, 10)

%spy(A)
eig(A)

V = chol(A)






%5.4
A = [1 1e10 1 1; 1e10 1 1 1e10; 1 1 1e-10 1; 1 1e10 1 1e10]
b=[1e10+3 2*1e10+2 3 1e-10 2*1e10+2]'
x=[1 1 1 1]'

b= A*x
[L, U, P]= lu(A)
y=L\(P*b)
x1=U\y

eig(A)


norm(b-A*x, inf)/norm(b,inf) %small residual 
cond(A, inf)


