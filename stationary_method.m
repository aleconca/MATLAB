function [x, iter, incr] = stationary_method(B, g, x0, tol, max_it)
% stationary_method Solve a linear system with a stationary iterative method of the form
% x^(k+1) = B*x^(k) + g
%
%   [x, iter, incr] = stationary_method(B, g, x0, tol, max_it)
%
%   Inputs  : B      = iteration matrix
%             g      = iteration vector
%             x0     = starting point
%             tol    = requested tolerance
%             max_it = maximum number of allowed iterations
%   Outputs : x      = approximation of the solution
%             iter   = perfomed iterations
%             incr   = last increment (can be used to estimate the error)

% Initialization
iter = 0;
incr = inf;
x    = x0;

while (incr > tol) && (iter < max_it)
  
  iter = iter+1;
  x_old = x;                   % Save the previous approximation

  x = B*x + g; %xnew = B*xold + g

  incr = norm(x - x_old, 2);   % 2-norm of the increment, you can choose the residual (as indicated in the text)
  % res = norm(b - A*x, 2);
end

end