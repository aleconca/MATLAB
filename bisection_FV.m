function [x,x_iter,n] = bisection_FV(f,a,b,tol)

    if f(a)*f(b) < 0
        
        Nmax = ceil(log((b-a)/tol)/log(2));
        n = 0;
        while n < Nmax
            n = n + 1;
            x_iter(n) = (a+b)/2;
            if f(x_iter(n)) == 0
                x = x_iter(n);     
            elseif f(a)*f(x_iter(n)) < 0
                b = x_iter(n);
            else
                a = x_iter(n);
            end

        end
        x = x_iter(n);
        
    end
    
end