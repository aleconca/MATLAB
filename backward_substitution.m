function x = backward_substitution(A,b)
  %BACKWARD_SUBSTITUTION Solve an upper triangular system using backward 
  %substitution method.
  %   x = backward_substitution(A,b)
  %
  %   Inputs  : A = system coefficient matrix
  %             b = right-hand side vector
  %   Outputs : x = solution vector

  n = length(b);
  x = 0*b;
  x(n) = b(n) / A(n,n);
  
  for (k = n-1:-1:1)
    x(k) = (b(k)-A(k,k+1:n)*x(k+1:n)) / A(k,k);
  end

end
