function [x,x_iter]=bisection(f,a,b,tol)

Nmax = ceil(log((b-a)/tol)/log(2));

for i=1:Nmax
    x_iter(i)=(b+a)/2;
    if f(x_iter(i))*f(a)<0
        b=x_iter(i);
    else
        a=x_iter(i);
    end
end
x=x_iter(end);
end