%22-07-21

%
f =@(x) x.^3 + x.^2 -2*x+3;
a=-1;
b=2;
m=20;

51/4

composite_midpoint(f,a,b,m)
composite_trapezoidal(f,a,b,m)
composite_simpson(f,a,b,m)%best result-->exact result


%
f = @(x) (1-x).^(-3) .*cos(2*pi*x);
a=-1;
b=0;

%i) Provide the plot of the function f in I.
figure;
x=linspace(a,b,1000);
plot(x, f(x))
xlabel('x')
ylabel('f(x)')
title('Plot of f(x) = (1-x).^(-3) .* cos(2*pi*x)')
grid on
%ii) Compute the Lagrange interpolant polynomial based on 6 equally spaced nodes.
x_nodes=linspace(a,b,6);
y_nodes=f(x_nodes);
p=polyfit(x_nodes,y_nodes,6);

x_plot=linspace(a,b,1000);
%y_plot
interp=polyval(p,x_plot);

%iii) Compute the associated interpolation error by using the infinite norm. Is the error unifomly distributed in I?
%err=norm(f(x_plot)-interp,'inf')-->max

error = abs(f(x_plot) - interp);

%If the error in the infinity norm is uniform over the interval I, 
%then the maximum error should remain relatively constant across the interval.

figure;
plot(x_plot, error, 'r-', 'LineWidth', 1.5); % Plot the error
hold on;
plot(x_plot, interp, 'b--', 'LineWidth', 1.5); % Plot the Lagrange polynomial
xlabel('x');
ylabel('Value');
legend('Error', 'Lagrange Polynomial');
title('Error and Lagrange Polynomial');
grid on;




%
f =@(x) exp(-x + 2)-3-x.^2;
a=0;
b=1.5;
phi=@(x) -log(x.^2+3)+2;
nmax=100;
tol=1e-6;

alpha = 0.7356;%exact zero

%i)approximate the zero using fixed-point
[x,iter]= fixed_point(phi, alpha,tol,nmax);
numel(iter)

%ii)are the hyp of Ostrowski theorem satisfied? Yes
dphi=@(x) -2.*x/(x.^2+3);
abs(phi(alpha))<1





%
n = 53;
A=3*eye(n)-diag(ones(n-1,1),1)-0.5*diag(ones(n-1,1),-1)-0.1*diag(ones(n-5,1),5)-0.05*diag(ones(n-5,1),-5);
xex=ones(1,n);
b=A*xex';

%We approximate the solution to the system by using a stationary Richardson
%method. Now:

%i) Verify if the Jacobi method is convergent.
D=diag(diag(A));
E=-tril(A,-1);

Bj=D\(D-A);
gj=D\b;
max(abs(eig(Bj)))<1

%ii) Verify if the Gauss-Seidel method is convergent.
Bgs=(D-E)\(D-E-A);
ggs=(D-E)\b;
max(abs(eig(Bgs)))<1%ok

%iii) Determine the number of iterations required by the Jacobi method to converge after 
% setting the tolerance tol=1e-8, the maximum number of iterations nmax=200 and the initial guess x0=zeros(n,1).
tol=1e-8;
nmax=200;
x0=zeros(n,1);

[x,iter]=stationary_method(Bj,gj,x0,tol,nmax);
numel(iter)
%iv) Determine the number of iterations required by the Gauss-Seidel method to converge after setting the tolerance tol=1e-8,  
% the maximum number of iterations nmax=200 and the initial guess x0=zeros(n,1).
tol=1e-8;
nmax=200;
x0=zeros(n,1);

[x,iter]=stationary_method(Bgs,ggs,x0,tol,nmax);
numel(iter)

