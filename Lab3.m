
%3.1 fixed point 

clear all
close all
clc

xi=sqrt(5)

% Method 1.

phi1  = @(x) 5 + x - x.^2;
dphi1 = @(x) 1 - 2*x;
abs(dphi1(xi))

% The absolute value of the derivative of phi at xi
% is greater than 1: the method will not converge.


% Method 2.

phi2  = @(x) 5./x;
dphi2 = @(x) -5./x.^2;
abs(dphi2(xi))

% The absolute value of the derivative of phi at xi
% is equal to 1: no theoretical conclusion can be stated in this case.


%4
phi4  = @(x) 1/2*(x + 5./x);
dphi4 = @(x) 1/2 - 5/2*(1./x.^2);
abs(dphi4(xi))

%the first derivative is = 0 , < 1 so it will converge 
% so we calculate the second one which will be != 0-> convergence rate is 2

%now we apply the fixed point method

tol= 1e-6;
maxit=1000;

x0=xi+0.001; % un intorno del punto fisso
[xi1, x_iter1] = fixed_point(phi1, x0, tol, maxit);
xi1
iter1=numel(x_iter1)-1

x0=xi+0.001; % un intorno del punto fisso
[xi2, x_iter2] = fixed_point(phi2, x0, tol, maxit);
xi2
iter2=numel(x_iter2)-1


x0=xi+0.001; % un intorno del punto fisso
[xi4, x_iter4] = fixed_point(phi4, x0, tol, maxit);
xi4
iter4=numel(x_iter4)





