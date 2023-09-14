function [x,x_iter]=newton(f,df,x0,tol,Nmax)

i=1;
err=1+tol;
% or we can set
% err = inf;
x_iter(i)=x0;

while i<=Nmax && err>tol
    
    if(abs(df(x_iter(i))) < 1e-8) % check if the derivative (denominator) is zero (or very small)
        break;
    end
        
    x_iter(i+1)=x_iter(i)-f(x_iter(i))/df(x_iter(i));
    err=abs(x_iter(i+1)-x_iter(i));
    i=i+1;
end

x=x_iter(end);
end