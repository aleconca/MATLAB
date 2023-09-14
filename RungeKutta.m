
clc
clear all
close all

% Set the Runge function
f = @(x) 1./(1 + x.^2); 
% Set the point to build the plot of the Rung function
xx = linspace(-5, 5, 1000);
fxx = f(xx);

% Plot the Runge function
figure
subplot(1,2,1)
hold on, box on
plot(xx, fxx, 'k-', 'LineWidth',2)
axis([-5.1 5.1 -0.4 1.2])
set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)
xlabel('x','FontSize',16)
ylabel('f(x)','FontSize',16)

% Plot the Linear interpolant polynomial
for k = [2:2:10]
  x = linspace(-5, 5, k);
  fx = f(x);
  coef = polyfit(x, fx, k-1);
  yy = polyval(coef, xx);
  plot(xx, yy, 'r-', 'LineWidth',2)
  norm((fxx-yy),'inf')
  pause
end

% Plot the Runge function
subplot(1,2,2)%impila grafici
hold on, box on%aggiungi plot a figura esistente creando un riquadro
plot(xx, fxx, 'k-', 'LineWidth',2)
axis([-5.1 5.1 -1.2 2.4])
set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)
xlabel('x','FontSize',16)
ylabel('f(x)','FontSize',16)

% Plot the Linear interpolant polynomial
for k = [12:2:16]
  x = linspace(-5, 5, k);
  fx = f(x);
  coef = polyfit(x, fx, k-1);    
  yy = polyval(coef, xx);
  plot(xx, yy, 'r-', 'LineWidth',2)
  norm((fxx-yy),'inf')
  pause
end

% The interpolation polynomials approach the function in the middle
% of the interval, but close to the boundaries increasing oscillations appear.

