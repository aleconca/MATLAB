% Authors A.G.Mauri

% This script computes the integral I(f) = \int_a^b f(x) dx 
% using the composite Gauss-Legendre quadrature formula

% Integration limits defining the integration interval
a = 0;
b = pi;

% uniform partition into K elements (users can modify the the number
% that determine the partitioning)  
K = 2;
% mesh size
h = (b-a)/K;
% partition of [a,b]
x = [a:h:b];

% choice of the number of the point in the Gauss formula
n = 5;

% Tabular of nodes and weights on the reference interval [-1, 1]
if (n==1)
   xi = 0;
   w  = 2;
elseif (n==2)
   xi = [-1/sqrt(3), 1/sqrt(3)];
   w  = [1, 1];
elseif (n==3)
   xi = [-sqrt(3/5), 0, sqrt(3/5)];
   w  =[5/9, 8/9, 5/9];
elseif (n==4)
   xi = [-sqrt(3/7+2/7*sqrt(6/5)), -sqrt(3/7-2/7*sqrt(6/5)), sqrt(3/7-2/7*sqrt(6/5)), sqrt(3/7+2/7*sqrt(6/5))];
   w  =[(18-sqrt(30))/36, (18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36];
elseif (n==5)
   xi = [-1/3*sqrt(5+2*sqrt(10/7)),-1/3*sqrt(5-2*sqrt(10/7)), 0, 1/3*sqrt(5-2*sqrt(10/7)), 1/3*sqrt(5+2*sqrt(10/7))];
   w  =[(322-13*sqrt(70))/900, (322+13*sqrt(70))/900, 128/225, (322+13*sqrt(70))/900, (322-13*sqrt(70))/900];
end

% initialize composite integral
IG = 0;
% composite  numerical integration
for i=1:K
% local change of coordinates from [-1,1] to [x_{i-1}, x_i]
    xa = x(i);
	xb = x(i+1);
% calculation of the quadrature nodes on I_i
    xn = h/2*xi + (xa+xb)/2;
    fn = sin(xn);
% calculation of the local weights
	wn = h/2*w;		 
    IG = IG + sum(wn.*fn);
end

IG

% Manual calculation of the exact integration of the given function in the
% integration interval
Iex = 2 ;

% error evaluation
err = abs(Iex - IG)