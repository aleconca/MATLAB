% Authors A.G.Mauri
% This script calculate and visualize the Lagrange Polynomial

close all
clear all

x0 = -1 ;
xn = 1 ;
L  = xn-x0 ;

% Set the number of the interpolation interval (point are equal n+1)
n  = 4 ;

% interpolation nodes
xi =[x0:L/n:xn];
   
% graphical representation to emulate continuous function
x_plt = [x0:L/1000:xn]';

% set number of the lagrange polynomial
% change the number from 0 to n with number<=n
number = 2 ;

% Plot the lagrange polynomial
for it=0:number
plot(x_plt,lagrange_polynomial_fun(xi,it+1,x_plt))
hold on
end
grid on

legend('\phi_0','\phi_1','\phi_2','location','south','orientation','horizontal')
print('ar_lagr_pol','-dpng')