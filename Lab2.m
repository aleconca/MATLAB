%Bisection method
%2.1
close all 
clear all
clc

f= @(x) x.^3-(2 + exp(1))*x.^2 + (2*exp(1) + 1)*x + (1 - exp(1)) - cosh(x - 1);

a=0.5;
b=5.5;

x_plot=linspace(a,b,1000);
plot(x_plot,f(x_plot));
grid on
pause

% from the plot two different roots can be identified. The bisection method, 
% nevertheless, can be applied only to the root on the right: the one
% belonging to the interval [3,5]. Actually, in the other case, there isn't
% any interval such that the first derivative has constant sign and the
% root belongs to this interval. Therefore the Bolzano theorem hypothesis
% are not respected in this case.

a=3;
b=5;
tol=1.e-3;
[x, x_iter]=bisection(f,a,b,tol);%beware, this computes the number of iterations using Nmax as a limit
x

% plot the iterations on the function

x_plot=linspace(a,b,1000);
figure
plot(x_plot,f(x_plot));
title("Iterations required to reach the tollerance:", length(x_iter))
grid on
hold on
disp('Press any key...');
for i=1:length(x_iter)
plot(x_iter(i),f(x_iter(i)),'r*')
pause
end

%approximate the roots
[x,x_iter,n] = bisection_FV(f,a,b,tol);








%2.2 Newton
close all 
clear all
clc

g= @(x) sin(x).*(1-x).^2;
dg=@(x) cos(x).*(1-x).^2-2*sin(x).*(1-x);
a=-0.5;
b=1.5;

x_plot=linspace(a,b,1000);
plot(x_plot,g(x_plot));
grid on
pause

%initial guess
x0=0.3;
tol=1.e-6;
Nmax=100;



[x1,x1_iter]=newton(g,dg,x0,tol,Nmax);
err1=abs(x1_iter-0)
x1

% estimate of the convergence order
p=log(err1(3:end)./err1(2:end-1))./log(err1(2:end-1)./err1(1:end-2));
figure;
plot(p);
title('convergence order, root x1=0');
% Simple zero => p = 2
% Multiple zero => p = 1=> Modified newton

%[x,x_iter] = modified_newton(g,dg,x0,tol,Nmax,m)

% Errors behaviour
figure;
semilogy(err1);
hold on
%semilogy(err2,'g');
%legend('x1=0','x2=1');
%title('Errors');

