%
f=@(x) cos(x);
a=0;
b=pi;
n = 11;

x_nodes=linspace(a,b,n);
y_nodes=f(x_nodes);
p=polyfit(x_nodes,y_nodes,n-1);

x_plot=linspace(a,b,1000);
%y_plot
interp=polyval(p,x_plot);

figure;
plot(x_nodes, y_nodes, x_plot, interp)
xlabel('x');
ylabel('f(x)');

norm(f(x_plot)-interp,'inf')


%the error is very small meaning that we are obtaining a good apporoximation;



%
n=50;
A=eye(n) + 0.2*diag(ones(n-1,1),+1) + 0.1*diag(ones(n-1,1),-1);

%ok N.S. condition
for i=1:n-1
    det(A(1:i,1:i))~=0;
end


val=0;

%A(i,j)
for j=1:n
    for i=1:n
        val=val+A(i,j); 
    end
    abs(A(j,j))>=val;
end

for i=1:n
    for j=1:n
        val=val+A(i,j); 
    end
    abs(A(i,i))>=val;
end

xex = [1:50]';
b=A*xex;
[L,U]=lu(A);

y=forward_substitution(L,b);
x=backward_substitution(U,y);
x

norm(xex-x)

%The presence of this small error can be attributed to round-off errors and
% numerical instability in the computations involved in LU factorization. 
% These errors can accumulate and affect the final result, even though the original problem is well-posed.

%Round-off errors occur due to limited precision in representing real numbers on a computer. 
% Operations such as addition, subtraction, and multiplication can introduce small errors due 
% to the truncation or rounding of digits. These errors can propagate through the LU factorization
% process and result in a small discrepancy between the computed solution and the exact solution.

%Numerical instability can arise when the matrix A is ill-conditioned or nearly singular. 
% Ill-conditioned matrices can amplify the effect of round-off errors and make the computations less accurate. 
% This can lead to a small error in the solution even when the exact solution is known.





