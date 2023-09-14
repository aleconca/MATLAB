%3.1
%a
phi = @(x) x-x.^3;
xi=0; %fixed point

%dphi = @(x) 1-3*(x.^2); %this is = 1 at zero so we cannot conclude anything

x0 = -1 + 2*rand(1)

tol=1e-6;
maxit=10000; %it coonverges but very slowly

[xi1, x_iter1] = fixed_point(phi, x0, tol, maxit);
xi1
iter = numel(x_iter1)-1

%2
phi = @(x) x+x.^3;
xi=0; %fixed point

x0 = -1 + 2*rand(1)

tol=1e-6;
maxit=10000; 

[xi1, x_iter1] = fixed_point(phi, x0, tol, maxit);
xi1
iter = numel(x_iter1)-1

%this diverges

%3.2
f=@(x) exp(x.^2).*log(x+1)-1;
tol=1e-3;
maxit=1000;

rootfinding_function_plot(f,-1,2);




