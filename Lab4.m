%4.1

%a
A = [1 2 3 4; 0 1 2 3; 0 0 1 2; 0 0 0 1];
b = [1 1 1 1]';

x = backward_substitution(A,b);
x

A\b

%b
At = [1 2 3 4; 0 1 2 3; 0 0 1 2; 0 0 0 1]';
b = [1 1 1 1]';

x = forward_substitution(A,b);
x

%c

At*A

y = forward_substitution(At,b)
x = backward_substitution(A,y)

(At*A)\b


%4.2

A = [2 10 4 0; 
     1 0 2 2; 
     1 4 0 2; 
     1 2 1 1];
b = [10 1 3 3]';

%[L, U, P]=lu(A)
[L, U]=lu(A)


y = forward_substitution(L,b)
x = backward_substitution(U,y)
A\b

%A non singular
n=4
det(A(1:n-2,1:n-2))
det(A(1:n-1,1:n-1))

%A is non singular, thus det(A)=det(L)*det(U)

detL=prod(diag(L))
detU=prod(diag(U))
detA=detL*detU




%4.3

clear all
close all
clc

A = [50 1 3; 1 6 0; 3 0 1];

[L, U] = lu(A);

x=[1 1 1]
b=x*A

%y = forward_substitution(L,b)
%x = backward_substitution(U,y)

%x

%A\b


%4.4

[L,U]=lu(A);
n=length(A);
B= A;


for(i=1 : 1 : n)
    I = zeros(n,1);
    I(i)=1
    y=forward_substitution(L,I)
    x=backward_substitution(U,y)

    B(1:n,i) = x
end


B

inv(A)




