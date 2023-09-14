function [xi, x_iter] = newton25(f, df, x0, tol, maxit, multiplicity)
  %NEWTON Find a root of the equation f(x) = 0 using the Newton method, starting from the initial guess x0.
  %
  %   [xi, x_iter] =  NEWTON(f, df, x0, tol, maxit, multiplicity)
  %
  %   Inputs  : f     = function handle to the function f(x)
  %             df    = function handle to the derivative of the function f(x)
  %             x0    = initial guess
  %             tol   = requested tolerance
  %             maxit = maximum number of iterations
  %             multiplicity = multiplicity of the root
  %   Output  : 
  %			xi = approximation of the root
  %			x_iter     = vector of the approximations of the root at each step

  if (nargin < 6)
	multiplicity = 1;
  end
  
  x_iter(1) = x0;

  for (iter = 1:maxit)
    newton_method = @(x) x - multiplicity*f(x)/df(x);
    x_iter(iter+1) = newton_method(x_iter(iter));

    if (abs (x_iter(iter+1) - x_iter(iter)) < tol)
      break;
    end

  end
  
  xi = x_iter(end);

end