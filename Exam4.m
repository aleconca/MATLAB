%15-07-21

%
n=75;
a=1:n;
A=diag(a)-diag(ones(n-2,1),2);

b=ones(n,1);

[L,U,x] =thomas(A,b);
x

%
f=@(x) sin(x).*sqrt(x+0.5);
a=-0.5;
b=1;

x_nodes=linspace(a,b,6);
y_nodes=f(x_nodes);
p=polyfit(x_nodes,y_nodes,6);

x_plot=linspace(a,b,1000);
interp_plot=polyval(p,x_plot);

figure;
error=abs(f(x_plot)-interp_plot);
plot(x_plot,error,'r-')
hold on;
plot(x_plot,interp_plot,'b--')
grid on;


%
Nmax=20;
tol=1e-6;
df= @(x) cos(x).*sqrt(x+0.5)+sin(x).*0.5.*sqrt(x+0.5)^(-1);
alpha=0;


%first check if one of the endpoints a or b is already equal to the zero. 
% If so, set that endpoint as the initial guess x0 for the Newton method. 
% Otherwise, it proceeds with the bisection method to approximate the zero, 
% and then uses the result as the initial guess for the Newton method.
[x,x_iter]=bisection(f,a,b,tol);
x
x_iter
%In the case of the function f(x) = sin(x)*sqrt(x+0.5) in the interval I = [-0.5, 1], 
% since there is a sign change within the interval, applying the bisection method will converge to one of the two zeros 
% within that interval. Which zero it converges to depends on the initial interval boundaries and the initial guess provided.

%[x,x_iter]=newton(f,df,x0,tol,Nmax);
%x
%To approximate multiple zeros of a function, you would need to apply the Newton's method with different initial guesses, 
% each one close to a different zero. 
% This way, you can iteratively approximate each zero separately.
%In summary, if a function has multiple zeros, the Newton's method applied with a specific initial guess 
% will only give you an approximation of one of those zeros, not all of them simultaneously.


p=Stimap(x_iter,numel(x_iter));%convergence order



%
f=@(x) (x.^5-32).*log(x);
a=0.5;
b=3;
m=20;%subintervals

composite_midpoint(f,a,b,m)%
composite_trapezoidal(f,a,b,m)%
%The midpoint rule has an error estimate of O(h^2), where h is the width of each subinterval. 
% This means that as you decrease the width of the subintervals (increase the value of N), the 
% approximation should converge to the true value at a rate of O(h^2). So, if you compare the 
% approximations obtained with different values of N, you can expect the approximation to improve as N increases.

%The trapezoidal rule has an error estimate of O(h^2) as well. 
% Similarly, increasing the number of subintervals should lead to a more accurate approximation.

composite_simpson(f,a,b,m)
%Simpson's rule has an error estimate of O(h^4). 
% This means that it has the potential to provide more accurate results compared to the midpoint and trapezoidal rules, 
% especially when the function being integrated is smooth and can be well approximated by a quadratic polynomial over each 
% pair of subintervals. Again, increasing the number of subintervals should improve the approximation.





%
n=50;
a=1:n;
A=diag(a)-diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);
xex=ones(n,1);
Nmax=200;
tol=1e-6;

D=diag(diag(A));
E=-tril(A,-1);
B=(D-E)\(D-E-A);

max(abs(eig(B)))<1;%convergent

b=A*xex;
g=(D-E)\b;

x0 = zeros(n, 1); 

[x, iter, incr] = stationary_method(B, g, x0, tol, Nmax);
numel(iter)

%By solving the system using the Gauss-Seidel method, you can verify if the method converges to the exact solution 
% within the specified tolerance. The solution error, computed as the difference between the obtained solution and 
% the exact solution, can give you an indication of the accuracy of the iterative method in this particular case.

% Check the solution error
solution_error = norm(x - xex, inf)%there is some discrepancy between the computed solution and the true solution.

%The error in the solution can arise due to several factors, such as: 
% the numerical precision of the iterative method, 
% the chosen tolerance, 
% the number of iterations performed, 
% or the specific properties of the system matrix A. 

% Even if you provide the exact solution xex to compute b, the iterative method may not converge exactly 
% to the true solution due to these factors.




