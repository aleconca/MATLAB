%LU with pivoting
%5.1

E= [4 1 1 1 5; 
    4 1 2 0 0;
    1 0 15 5 1;
    0 2 4 10 2;
    3 1 2 4 20];

b = [12 19 22 18 30]';

%1
%A non-singular-->NO, LU fact without Pivoting doesn't exist
det(E(1:length(E)-4,1:length(E)-4))
det(E(1:length(E)-3,1:length(E)-3))%det=0
det(E(1:length(E)-2,1:length(E)-2))
det(E(1:length(E)-1,1:length(E)-1))
% sum by row
for i=1:1:length(E)
   x=E(i,i)>sum(E(i,i:length(E)));
   x
end
% sum by columns
for j=1:1:length(E)
   x=E(j,j)>sum(E(j:length(E),j));
   x
end
% spd symmetric and all eig are positive
%E is not symmetric
eig(E)

%2
[L, U, P]=lu(E);
%spy(L)

y=forward_substitution(L,P*b)
x=backward_substitution(U,y)

y=L\(P*b)
x=U\y

E\b



%5.2 Cholesky
%A spd

A= [44 15 29 26 119;
    15 33 32 18 15;
    29 32 252 112 73;
    26 18 112 124 90;
    119 15 73 90 430];
%1
%A is symmetric
eig(A)

%2
%3
R=chol(A)
MyChol(A)'

%4
b = [1 1 1 1 1]';

y=R'\b
x=R\y

A\b



%5.3 Thomas

a=ones(10,1);
b=ones(10-1,1);
c=ones(10-1,1);

for i=1:1:10
    a(i)=i;
    if(i<10)
     b(i)=(11+i-1);
     c(i)=(102+i-1);
    end 
end


A = diag(a) + diag(b,1) + diag(c,-1)

S = sparse(A)

[Bout,id] = spdiags(A)

x=ones(10,1)
b=A*x
[L, U, x] = thomas(A,b)

A\b
