function [x, iter, incr] = prec_rich_method(A, b, P, alpha, x0, tol, max_it)
% prec_iter_method Solve a linear system using the preconditioned Richardson
% method of the form
% x^(k+1) = x^(k) + alpha * inv(P) * r^(k)
%
%   [x, iter, incr] = prec_rich_method(A, b, P, alpha, x0, tol, max_it)
%
%   Inputs  : A      = coefficient matrix
%             b      = right hand side vector
%             P      = preconditioner
%             alpha  = relaxation parameter
%             x0     = initial guess
%             tol    = requested tolerance
%             max_it = maximum number of allowed iterations
%   Outputs : x      = approximation of the solution
%             iter   = perfomed iterations
%             incr   = last increment (can be used to estimate the error)

% Initialization
iter = 0;
incr = inf;
x    = x0;
r    = b - A*x;

while (incr > tol) && (iter < max_it)
  
  iter = iter+1;
  
  z = P \ r;
  x = x + alpha * z; 
  r = r - alpha * A * z; % b - A*x
  
  incr = norm(alpha * z, 2);       % 2-norm of the increment x - xold
    
end

end