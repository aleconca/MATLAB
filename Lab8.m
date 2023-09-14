%8.1

f=@(x) x.*sin(x);
a=-2;
b=6;

%a)
%equally spaced nodes
%degree 4,6,8
n=4;
x_nodes=linspace(a,b,n+1);
y_nodes=f(x_nodes);
p=polyfit(x_nodes,y_nodes,n);

%plot the function
x_plot=linspace(a,b,1000);
y_plot=f(x_plot);
interp_plot=polyval(p,x_plot);

%figure;
plot(x_plot, f(x_plot), '--', x_plot, interp_plot, x_nodes, y_nodes, 'o', 'LineWidth', 2, 'MarkerSize', 8)





%b)
% calculation of the maximum in the derivative term 
if (n == 4)
  % the n+1-th derivative is 5 sin(x) + x cos(x). An upper bound of its absolute value is
  deriv_upper_bound = 11 ; %abs(f^n+1) is max at
elseif (n == 6)
  % the n+1-th derivative is -7 sin(x) - x cos(x). An upper bound of its absolute value is
  deriv_upper_bound = 13 ;
elseif (n == 8)
  % the n+1-th derivative is 9 sin(x) + x cos(x). An upper bound of its absolute value is
  deriv_upper_bound = 15 ;
else
  error('Wrong n');
end
% evaluation of the upper bound of the error 
err_est = (x_nodes(2)-x_nodes(1))^(n+1)/(4*(n+1)) * deriv_upper_bound;

% Calculation of the error
% err_norm = norm((f_plot - interp_plot),'inf');
err = abs(y_plot - interp_plot);

% Plotting the 
%figure;
plot(x_plot, err, '-b','LineWidth',2,'MarkerSize',12);
hold on;
plot(x_plot, err_est*ones(size(x_plot)), '-r','LineWidth',2,'MarkerSize',12)





%c)
n = 5; % n=10, n=20
xint = linspace(a, b, n+1);
yint = f(xint);
y_plot = interp1(xint, yint, x_plot);

%figure;
plot(x_plot, f(x_plot), 'k-', x_plot, y_plot, 'r-', xint, yint, 'rx', 'LineWidth', 2, 'MarkerSize', 8)
title('Piecewise Linear interpolation') ;

% c-2) Computation of the spine line interpolant using spline Matlab command  
% on equally spaced nodes
xsp = linspace(a, b, n+1);
ysp = f(xsp);
y_plot = spline(xsp, ysp, x_plot);

%figure;
plot(x_plot, f(x_plot), 'k-', x_plot, y_plot, 'r-', xsp, ysp, 'rx', 'LineWidth', 2, 'MarkerSize', 8)
title ('Spline interpolation')











%8.2
%a)see Runge Kutta.m


% b) Compute the piecewise linear interpolation (interp1 command) with
% n = 1, 2, 4, 8, 16, 32 elements and plot the result. Find the error

clear all
a = -5;
b = 5;
f = @(x) 1./(1+x.^2); 
x_plot = linspace(a, b, 101);

n_vect = [1 2 4 8 16 32 64 128];

figure
for (i = 1:numel(n_vect))
  n = n_vect(i)
  
  x = linspace(a, b, n+1);
  y = f(x);

  y_plot = interp1(x, y, x_plot);
  plot(x_plot, f(x_plot), 'k-', x_plot, y_plot, 'r-', x, y, 'rx', 'LineWidth', 2, 'MarkerSize', 8)
  axis([a-0.2 b+0.2 -0.1+min(f(x_plot)) 0.1+max(f(x_plot))])
  set(gca,'FontSize', 16)
  set(gca,'LineWidth', 1.5)
  
  error(i) = max(abs(f(x_plot) - y_plot));
  pause
end

H = 10./n_vect
pause

figure
loglog(H, error, 'bx-', 'LineWidth', 2, 'MarkerSize', 8)
hold on, box on
loglog(H, H.^(2), 'g-', 'LineWidth', 2)
axis([1e-1 1e1 1e-3 1e2])
set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)
xlabel('h','FontSize',16)
ylabel('error','FontSize',16)
legend('error', 'O(h^2)', 'Location', 'ne');

% From the graph we can conclude that the convergence order is 2,
% as predicted by the theory.
% Calculation of the error
order = (log(error(1:end-1) ./ error(2:end))/log(2))'
% or using the diff command build in in matlab
p = -diff(log(error)) / log(2)

% c) Interpolation of the Runge function using a piecewise cubic spline on
% the same number of iterpolating points

figure
for (i = 1:numel(n_vect))
  n = n_vect(i)
  
  x = linspace(a, b, n+1);
  y = f(x);

  y_plot = spline(x, y, x_plot);
  plot(x_plot, f(x_plot), 'k-', x_plot, y_plot, 'r-', x, y, 'rx', 'LineWidth', 2, 'MarkerSize', 8)
  axis([a-0.2 b+0.2 -0.1+min(f(x_plot)) 0.1+max(f(x_plot))])
  set(gca,'FontSize', 16)
  set(gca,'LineWidth', 1.5)
  
  pause
end

% d) Computation of the linera interpolant by using the --->Chebyshev nodes

xx = linspace(-5, 5, 1000);

figure;
for (i = 1:numel(n_vect)-2)
  n = n_vect(i) ;
  
  % Calculation of the Chebyshev nodes
  ii = 0:n;
  x_cap = -cos(pi*ii/n);
  x = 0.5*(a+b) + 0.5*(b-a)*x_cap;
  y = f(x);
  
  coef = polyfit(x, y, n);    
  yy = polyval(coef, xx);
  
  plot(xx, f(xx), 'k-', xx, yy, 'r-', x, y, 'rx', 'LineWidth', 2, 'MarkerSize', 8)
  axis([a-0.2 b+0.2 -0.1+min(f(xx)) 0.1+max(f(xx))])
  set(gca,'FontSize', 16)
  set(gca,'LineWidth', 1.5)
  
  error(i) = max(abs(f(xx) - yy));
  
  
  pause
  
end







%8.3

%a)
clc
clear all
close all

% Setting the function and its visualization
a = -1;
b = 1;
f = @(x) abs(x - pi/12); 

xx = linspace(a, b, 1000);
fxx = f(xx);

figure
subplot(1,2,1)
hold on, box on
plot(xx, fxx, 'k-', 'LineWidth',2)
set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)
xlabel('x','FontSize',16)
ylabel('f(x)','FontSize',16)

% a) The interpolating polynomials Π n f (x ) using equally
% spaced nodes show the same ---->Runge-s phenomenon

for k = [2:2:10]
  x = linspace(a, b, k);
  fx = f(x);
  coef = polyfit(x, fx, k-1);
  yy = polyval(coef, xx);
  plot(xx, yy, 'r-', 'LineWidth',2)
  axis([a-0.2 b+0.2 -0.1+min(f(xx)) 0.1+max(f(xx))])
  pause
end

subplot(1,2,2)
hold on, box on
plot(xx, fxx, 'k-', 'LineWidth',2)
set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)
xlabel('x','FontSize',16)
ylabel('f(x)','FontSize',16)

for k = [12:2:16]
  x = linspace(a, b, k);
  fx = f(x);
  coef = polyfit(x, fx, k-1);    
  yy = polyval(coef, xx);
  plot(xx, yy, 'r-', 'LineWidth',2)
  axis([a-0.2 b+0.2 -0.1+min(f(xx)) 0.1+max(f(xx))])
  pause
end

% The interpolation polynomials approach the function in the middle of the 
% interval, but close to the boundaries increasing oscillations appears.



%d)least square
%Interpolation of the points generated with the above function with nodes n = 10
% and using the least square approach with varying the degree of the interpolant.

close all

a = -1;
b = 1;
f = @(x) abs(x - pi/12); 
x_plot = linspace(a, b, 101);
x_nodes = linspace(a, b, 10);
y_nodes = f(x_nodes);

plot(x_nodes, y_nodes, 'rO');
hold on;

for ii = 2:1:15
    v = polyfit(x_nodes, y_nodes, ii);
    y_plot = polyval(v, x_plot);
    plot(x_plot, y_plot)
    pause
end
