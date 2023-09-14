%07-09-21
f=@(x) cos(x)+(x-pi./2).^2;
a=0;
b=pi/2;

n=5;

x_nodes=linspace(a,b,n+1);
y_nodes=f(x_nodes);
p=polyfit(x_nodes,y_nodes,n);

x_plot=linspace(a,b,1000);
%y_plot=f(x_plot);
interp=polyval(p,x_plot);


% Plot the function and the interpolating polynomial
figure
plot(x_plot, f(x_plot), 'b-', 'LineWidth', 2)
hold on
plot(x_plot, y_plot, 'r--', 'LineWidth', 2)
legend('f(x)', 'Interpolating Polynomial')
xlabel('x')
ylabel('f(x)')
title('Interpolating Polynomial of Degree 5')
grid on


%
n=20;
composite_midpoint(f,a,b,n)
composite_trapezoidal(f,a,b,n-1)
composite_simpson(f,a,b,(n-1)/2)


%
n=80; 
a=7*ones(n,1);
A=diag(a)-3*diag(ones(n-1,1),1)-3*diag(ones(n-1,1),-1);
b=A*ones(n,1);

[L,U,x]=thomas(A,b);
y=forward_substitution(L,b);
x1=backward_substitution(U,y);

x=A\b;

norm(x-x1)

%the norm of the error vector is approximately 4.9051e-15, which is a very small value close to zero, 
%it suggests that the computed solution is very accurate and close to the exact solution.


%
f=@(x) sin(x-0.5)*log(x+5)^2;
a=-4;
b=4;

for i=2:1:6
    comp_trap(i)= composite_trapezoidal(f,a,b,2^i)
end   

errT=abs(3.7518-comp_trap)
p=-diff(log(errT))/log(2)


%
f=@(x) sin(x)*sqrt(x+0.5);
a=-0.5;
b=1;

nmax=20;
tol=1e-06;
df=@(x) cos(x)*sqrt(x+0.5)+0.5*sin(x)*sqrt(x+0.5).^(-1);

x0=bisection(f,a,b,tol);
%newton(f,df,x0,tol,nmax) simple zero

for m=1:1:6
 c=modified_newton(f,df,x0,tol,nmax,m);
 p = stimap(approximation, f, x0)<=2
end

