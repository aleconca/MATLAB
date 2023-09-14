%
f=@(x) sin(x).*sqrt(x+0.5);
a=-0.5;
b=1;
tol=1e-6;
Nmax=20;
df=@(x) cos(x).*sqrt(x+0.5)+sin(x)*0.5*sqrt(x+0.5).^(-1);

x0=bisection(f,a,b,tol);
[x,x_iter]=newton(f,df,x0,tol,Nmax);
numel(x_iter)

err1=abs(x_iter-0)

% estimate of the convergence order
p=log(err1(3:end)./err1(2:end-1))./log(err1(2:end-1)./err1(1:end-2));
p




%
f =@(x) x.^2 .* (x - 3);
Nmax=200;
a=-1;
b=5;
df=2.*x.*(x-3)+x.^2;

x0=2;
%x1=newton(f,df,x0,tol,Nmax);
x0=0.5;
%x2=modified_newton(f,df,x0,tol,Nmax,2)



%
A=[30 5 30 5 5; 5 55 5 30 5; 30 5 55 5 30; 5 30 5 55 5; 5 5 30 5 30];
b=ones(5,1);
%A=A' symmetry
eig(A)>0  %positive-def

H=chol(A);%upper
H_T=H';%lower
y=forward_substitution(H_T,b);
x=backward_substitution(H,y);


%
n=75;
a=1:n;
A=diag(a)-diag(ones(n-2,1),2)

b=ones(n,1);

[L,U,x] =thomas(A,b);
x


%

n = 10;
A = 2*diag(ones(n, 1)) - 1*diag(ones(n-1, 1), -1) - 1*diag(ones(n-1, 1), 1);

b = ones(n, 1); 
%Compute the solution x in the following scenarios:

%i)  P = diag(diag(A)), and α chosen as the optimal accelaration parameter;
%ii) P = diag(diag(A)), α = 0.5;
%iii) P = eye(n, n), α = 0.4.

x0 = zeros(n, 1);
tol = 1e-6; 
maxit = 1000; 

%for all the three cases. Comment on the obtained results.

%i)
P1 = diag(diag(A));
alpha_opt=2/(max(abs(eig(inv(P)*A)))+min(abs(eig(inv(P)*A))));
[x1, iter, incr] =prec_rich_method(A,b,P1,alpha_opt,x0,tol,maxit);
x1
numel(iter)
incr


%ii)
P2 = diag(diag(A));
alpha2=0.5;
[x2, iter, incr] =prec_rich_method(A,b,P2,alpha2,x0,tol,maxit);
x2
numel(iter)
incr

%iii)
P3 = eye(n, n);
alpha3=0.4;
[x3, iter, incr] =prec_rich_method(A,b,P3,alpha3,x0,tol,maxit);
x3
numel(iter)
incr


%
%Consider the function f(x) = (1-x).^(-3) .*cos(2*pi*x) in the interval I=[-1 ,0].
f=@(x) (1-x).^(-3) .*cos(2*pi*x);
a=-1;
b=0;
%i) Provide the plot of the function f in I.
x_plot=linspace(a,b,1000);
y_plot=f(x_plot);

figure;
plot(x_plot,y_plot)
xlabel('x');
ylabel('f(x)');

%ii) Compute the Lagrange interpolant polynomial based on 6 equally spaced nodes.
x_nodes=linspace(a,b,6);
y_nodes=f(x_nodes);
p=polyfit(x_nodes,y_nodes,6);

interp_plot=polyval(p,x_plot);

%iii) Compute the associated interpolation error by using the infinite norm. Is the error unifomly distributed in I?
norm(y_plot-interp_plot,'inf')
error=abs(y_plot-interp_plot);

plot(x_plot, error, 'r-'); % Plot the error
hold on;
plot(x_plot, interp_plot, 'b--'); % Plot the Lagrange polynomial
grid on;



%
%Consider the function f(x)=sin(x)*sqrt(x+0.5) in the interval I=[-0.5,1] 
% and the associated Lagrange interpolating polynomial, L, based on 6 equispaced interpolation nodes. 

f=@(x) sin(x).*sqrt(x+0.5);
a=-0.5;
b=1;

x_nodes=linspace(a,b,6);
y_nodes=f(x_nodes);
p=polyfit(x_nodes,y_nodes,6);

x_plot=linspace(a,b,1000);
interp_plot=polyval(p,x_plot);

figure;
error=abs(f(x_plot)-interp_plot);
plot(x_plot,error,'r-')
hold on;
plot(x_plot,interp_plot,'b--')
grid on;

%Lagrange interpolant theoretical error:
%h=(b-a)./6;
%H=h.^7;
%df7=@(x)
%err= H/(4*7)*max(abs(df7(0.4)));
%err

%Chebshev nodes
f=@(x) 1./(1+x.^2);
n=[1 2 4 8 16 32];

for i=1:numel(n)

    niter=n(i);

    ii=0:niter;
    xcap=-cos((pi*ii)/niter);
    x=(a+b)./2 + (b-a)./2*xcap;
    y=f(x);

   
    p=polyfit(x,y,niter);

    x_plot=linspace(-5,5,1000);
    interp=polyval(p,x_plot);

    %figure;
    %plot(x_plot,interp)

    error(i)=max(abs(f(x_plot)-interp));
end

order= -diff(log(error))/log(2)




%least squares
f=@(x) abs(x-pi/12);
a=-1;
b=1;
n=10;

for m=2:2:6
   x_nodes=linspace(a,b,10);
   y_nodes=f(x_nodes);
   p=polyfit(x,y,m);

   x_plot=linspace(a,b,1000);
   interp=polyval(p,x_plot);

   %figure;
   %plot(x_plot,interp)
end




%
f =@(x) cos(x)+(x-pi/2).^2;
a=0;
b=pi/2;
m=20;

composite_midpoint(f,a,b,m)
composite_trapezoidal(f,a,b,m)
composite_simpson(f,a,b,m)


%
f=@(x) sin(x-0.5)*log(x+5).^2;
a=-4;
b=4;

for i=1:6
    res(i)=composite_trapezoidal(f,a,b,2^i);
end

error=abs(3.7518-res);
p=-diff(log(error))/log(2);
p


%
p =@(x) x.^3 + 5*x.^2 - x + 4;
a=0;
b=3;
m=15;

composite_midpoint(p,a,b,m)
composite_trapezoidal(p,a,b,m)
composite_simpson(p,a,b,m)



%
f=@(x) (x^5-32)*log(x);
a=0.5;
b=3;
m=20;

composite_midpoint(f,a,b,m)
composite_trapezoidal(f,a,b,m)
composite_simpson(f,a,b,m)




%
%Consider the linear system Ax = b, where A is the matrix defined as

n = 53;
A=3*eye(n)-diag(ones(n-1,1),1)-0.5*diag(ones(n-1,1),-1)-0.1*diag(ones(n-5,1),5)-0.05*diag(ones(n-5,1),-5);

%and b is chosen such that the exact solution is xex=ones(1,n).
xex=ones(1,n);
b=A*xex';

D=diag(diag(A));
E=-tril(A,-1);

%We approximate the solution to the system by using a stationary Richardson method. Now,
%i) Verify if the Jacobi method is convergent.
Bj=D\(D-A);
gj=D\b;

max(abs(eig(Bj)))<1

%ii) Verify if the Gauss-Seidel method is convergent.
Bgs=(D-E)\(D-E-A);
ggs=(D-E)\b;

max(abs(eig(Bgs)))<1

%iii) Determine the number of iterations required by the Jacobi method to converge after 
% setting the tolerance tol=1e-8, the maximum number of iterations nmax=200 and the initial guess x0=zeros(n,1).
tol=1e-8;
nmax=200;
x0=zeros(n,1);

[x, iter, incr] = stationary_method(Bj,gj,x0,tol,nmax);
numel(iter)

%iv) Determine the number of iterations required by the Gauss-Seidel method to converge 
% after setting the tolerance tol=1e-8,  the maximum number of iterations nmax=200 and the initial guess x0=zeros(n,1).
tol=1e-8;
nmax=200;
x0=zeros(n,1);

[x, iter, incr] = stationary_method(Bgs,ggs,x0,tol,nmax);
numel(iter)




