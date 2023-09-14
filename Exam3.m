%09-07-2020

%
f= @(x) cos(x).*sin(x);
a=0;
b=pi/2;
n=5;

x_nodes=linspace(a,b,n)
y_nodes=f(x_nodes)
p=polyfit(x_nodes,y_nodes,n)

x_plot=linspace(a,b,1000)
y_plot=f(x_plot)
interp_plot=polyval(p,x_plot)

figure;
plot(x_nodes, y_nodes,  x_plot, y_plot,  x_plot,interp_plot)


%
composite_trapezoidal(f,a,b,10)


%
A = 3*diag(ones(n, 1)) - 0.5*diag(ones(n-1, 1), -1) - 0.5*diag(ones(n-1, 1), 1)%tridiagonal,direct,thomas

%
f =@(x) 1/4 *x.^2 - x;
phi=@(x) 2*x - 1/4*x.^2;

alpha=4;
tol=1e-6;
nmax=200;
x0=0+rand(1);


[xi, x_iter] = fixed_point(phi,x0,tol,nmax)
numel(x_iter)
xi

%
f= @(x) cos(x).*sin(x);
a=0;
b=3/4*pi;
tol_b = 1e-2;
tol_n = 1e-6;

[x,x_iter]=bisection(f,a,b,tol_b);
df=@(x) cos(x).^2-sin(x).^2;
[xn,x_itern]=newton(f,df,x,tol_n,200)
numel(x_iter)
xn

%
p = @(x) x.^3 + 5*x.^2 - x + 4;
a=0;
b=3;
m=15;

composite_simpson(f,a,b,m)

%
f =@(x) x.^4 + 2*x.^2 - x - 3;
df=@(x) 4*x.^3 + 4*x -1;

phi1 =@(x) (3 + x - 2*x.^2).^(1/4);
phi2 =@(x) ((x + 3 - x.^4)/2).^(1/2);
phi3 =@(x) ((x+3)/(x.^2+2)).^(1/2);

alpha=1.1241;
tol=1e-6;
nmax=200;

[x1, x_iter1] = fixed_point(phi1,alpha,tol,nmax);%consistent
numel(x_iter1)
x1
[x2, x_iter2] = fixed_point(phi2,alpha,tol,nmax);
numel(x_iter2)
x2
[x3, x_iter3] = fixed_point(phi3,alpha,tol,nmax);%consistent
numel(x_iter3)
x3


%Richardson schemes
n = 10;
A = 2*diag(ones(n, 1)) - 1*diag(ones(n-1, 1), -1) - 1*diag(ones(n-1, 1), 1)
%eig(A)>0


b = ones(n, 1);

P = diag(diag(A));
%eig(P1)>0
alphaopt=2/(max(eig(inv(P)*A))+min(eig(inv(P)*A)));

alpha2 = 0.5;

P3 = eye(n, n); 
%eig(P3)>0
alpha3= 0.4;

x0 = zeros(n, 1);
tol = 1e-6; 
maxit = 1000;

%Comment on the obtained results.
[x1, iter1,incr1] = prec_rich_method(A, b, P, alphaopt, x0, tol, maxit);
iter1
[x2, iter2, incr2] = prec_rich_method(A, b, P, alpha2, x0, tol, maxit);
iter2
[x3, iter3, incr3] = prec_rich_method(A, b, P3, alpha3, x0, tol, maxit);
iter3


