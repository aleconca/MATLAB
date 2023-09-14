function [x,x_iter] = modified_newton(f,df,x0,tol,Nmax,m)

k = 0;
x_iter(k+1) = x0;
% or we can set
% err = 1 + tol;
err = inf;

while (err > tol) && (k <= Nmax)

    if(abs(df(x_iter(k))) < 1e-8) % check if the derivative (denominator) is zero (or very small)
        break;
    end
    
    k = k+1;
    x_iter(k+1) = x_iter(k) - m*f(x_iter(k))/df(x_iter(k));
    %err = abs( f(x_iter(k))/df(x_iter(k)) );
    err = abs(x_iter(k+1) - x_iter(k));
    
end
x = x_iter(end);


end