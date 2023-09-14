function [xi, x_iter] = fixed_point(phi, x0, tol, maxit)
  %FIXED_POINT Find a root of the equation f(x) = 0 using the fixed point
  %   method x = phi(x), starting from the initial guess x0.
  %
  %   [xi, x_iter] = FIXED_POINT(phi, x0, tol, maxit)
  %
  %   Inputs  : phi     = function handle to the iteration function
  %             x0    = initial guess
  %             tol   = requested tolerance
  %             maxit = maximum number of iterations
  %   Output  : 
  %			xi = approximation of the root
  %			x_iter     = vector of the approximations of the root at each step

  x_iter(1) = x0;

  for (iter = 1:maxit)

    x_iter(iter+1) = phi(x_iter(iter));

    if (abs (x_iter(iter+1) - x_iter(iter)) < tol)
      break;
    end

  end
  
  xi = x_iter(end);

end
