%8.1
%a)
f=@(x) exp(x);

a=-1;
b=1;


for n=[2:1:38]
x_nodes=linspace(a,b,n);
y_nodes=f(x_nodes);

x_plot=linspace(a,b,1000);
y_plot=f(x_plot);
p=interp1(x_nodes,y_nodes,x_plot);

max(abs(y_plot-p))<1e-3
n %n=38
end



%8.2
clc
clear all
close all

f = @(x) 1./((x - 0.3).^2 + 0.01) + 1./(x.^2 + 0.04) - 6;

% a)

n = 10;
x = linspace(-1, 3, n)';
fx = f(x);

coeff_polyfit = polyfit(x, fx, n-1);

x_plot = linspace(-1, 3, 1000);

p_equally = polyval(coeff_polyfit, x_plot);

figure
subplot(2,1,1)
plot(x_plot, f(x_plot), '-b', 'LineWidth', 2)
hold on, box on
plot(x_plot, p_equally, '-r', 'LineWidth',2)
plot(x, f(x), 'rx','LineWidth',2, 'MarkerSize',10)


% b)
err_equally = abs(f(x_plot) - p_equally);

subplot(2,1,2)
hold on, box on
plot(x_plot, err_equally, '-r', 'LineWidth',2)




%8.3
clc
clear all
close all

f = @(x) x.^5 + 1;

x = [0 0.2 0.4 0.6];
p3 = polyfit(x, f(x), 3)

x_plot = linspace(0, 0.6, 1000);
p3_plot = polyval(p3, x_plot);

figure
plot(x, f(x), 'rx','LineWidth', 2, 'MarkerSize', 10)
hold on, box on
plot(x_plot, f(x_plot), 'b--', 'LineWidth', 2)
plot(x_plot, p3_plot, 'r-', 'LineWidth', 2)
axis([-0.05 0.65 0.99 1.09])
set(gca, 'FontSize', 16)
set(gca, 'LineWidth', 1.5)
xlabel('x','FontSize',16)
ylabel('f(x)','FontSize',16)
leg = legend('nodes', 'f', 'p3', 'Location', 'nw');

err = max(abs(f(x_plot) - p3_plot))





%8.4
clc
clear all
close all

x_plot = linspace(0, 2, 1000);
f = @(x) 1*(x >= 1);

figure
plot(x_plot, f(x_plot), 'b-', 'LineWidth', 2)
hold on, box on
axis([-0.1 2.1 -1 2])
set(gca, 'FontSize', 16)
set(gca, 'LineWidth', 1.5)
xlabel('x','FontSize',16)
ylabel('f(x)','FontSize',16)

% Equispaced nodes

for n = [2, 4, 8, 16, 32]
  x = linspace(0, 2, n);

  coeff = polyfit(x, f(x), n-1);

  pn = polyval(coeff, x_plot);

  plot(x_plot, pn, 'r-', 'LineWidth',2)
  pause
end

% The result is very poor, for many reasons.
% The function to be interpolated is discontinuous,
% the interpolating polynomial has too many oscillations and
% the computation of the polynomial becomes ill-conditioned when
% the number of nodes increases.

% Piecewise linear polynomials

figure
plot(x_plot, f(x_plot), 'b-', 'LineWidth', 2)
hold on, box on
axis([-0.1 2.1 -1 2])
set(gca, 'FontSize', 16)
set(gca, 'LineWidth', 1.5)
xlabel('x','FontSize',16)
ylabel('f(x)','FontSize',16)

for n = [2, 4, 8, 16, 32]
  x = linspace(0, 2, n);

  s_plot = interp1(x,f(x),x_plot);
  plot(x_plot, s_plot, 'g-','LineWidth', 2)

  pause
end

% With linear splines the interpolation approximates the function
% with increasing precision, since for big n the slope at
% the discontinuity gets higher.

figure
plot(x_plot, f(x_plot), 'b-', 'LineWidth', 2)
hold on, box on
axis([-0.1 2.1 -1 2])
set(gca, 'FontSize', 16)
set(gca, 'LineWidth', 1.5)
xlabel('x','FontSize',16)
ylabel('f(x)','FontSize',16)

x = linspace(0, 2, 15);

for n = [2, 4, 8, 16, 32]
  
  v = polyfit(x,f(x), n);
  s_plot = polyval(v, x_plot);
  plot(x_plot, s_plot, 'g-','LineWidth', 2)

  pause
end

% Once again the result is very poor.
% The function to be interpolated is discontinuous,
% the interpolating polynomial has too many oscillations and
% the computation of the polynomial becomes ill-conditioned when
% the number of nodes increases because the interpolant tends to be smooth.

