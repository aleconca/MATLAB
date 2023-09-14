%choice for the alpha parameter in R.S.
%7.2

A=4*eye(50)-1*diag(ones(49,1),1)-1*diag(ones(49,1),-1)-1*diag(ones(48,1),2)-1*diag(ones(48,1),-2);
T=2*eye(50)-1*diag(ones(49,1),1)-1*diag(ones(49,1),-1);

%a ok spd
eig(A)>0
eig(T)>0

%b
x0=zeros(50,1);
b=0.2*ones(50,1);
tol=1e-6;
alpha1=0.2;
alpha2=0.33;
max_it=10000;
P=eye(50);

%inv(P)*A==A

%max(abs(eig(I - alpha1*A)))
%max(abs(eig(I - alpha2*A)))

[x, iter, incr] = prec_rich_method(A, b, P, alpha1, x0, tol, max_it)
[x, iter, incr] = prec_rich_method(A, b, P, alpha2, x0, tol, max_it)%no convergence

%alpha1 3355 iter

alphaopt= 2/(max(eig(inv(P)*A))+min(eig(inv(P)*A)))
[x, iter, incr] = prec_rich_method(A, b, P, alphaopt, x0, tol, max_it)

%c
P=T;
[x, iter, incr] = prec_rich_method(A, b, P, alphaopt, x0, tol, max_it)


cond(A)
cond(inv(P)*A) %way smaller ->
% Preconditioning the linear system with T results both in a
% better performance of the method (less iterations are needed) and a smaller condition number.

Tinv_A = T\A %inv(T)*A but more performant
