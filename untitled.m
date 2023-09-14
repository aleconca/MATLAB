clc
clear all
close all

x_plot = linspace(0, 2, 1000);
f = @(x) 1*(x >= 1);

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

