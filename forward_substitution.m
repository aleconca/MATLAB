function x = forward_substitution(A,b)
  %FORWARD_SUBSTITUTION Solve an lower triangular system using forward 
  %substitution method.
  %   x = forward_substitution(A,b)
  %
  %   Inputs  : A = system coefficient matrix
  %             b = right-hand side vector
  %   Outputs : x = solution vector

  n = length(b);
  % x = zeros(size(b));
  x = 0*b;
  x(1) = b(1) / A(1,1);
  
  for i=2:n
    x(i) = (b(i) - A(i, 1:i-1) * x(1:i-1)) / A(i,i);
  end

end
