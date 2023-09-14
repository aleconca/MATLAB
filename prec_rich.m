function [x, iter] = prec_rich(A, b, P, x0, alpha, tol, maxit)
    % Preconditioned Richardson method for solving Ax = b
    % Inputs:
    %   A: coefficient matrix
    %   b: right-hand side vector
    %   P: preconditioning matrix
    %   x0: initial guess
    %   alpha: acceleration parameter
    %   tol: tolerance for convergence
    %   maxit: maximum number of iterations
    % Outputs:
    %   x: solution vector
    %   iter: number of iterations taken
    
    x = x0;
    iter = 0;
    residual_norm = inf;
    
    while (residual_norm > tol) && (iter < maxit)
        iter = iter + 1;
        
        % Compute residual vector
        r = b - A * x;
        
        % Preconditioning step
        z = P \ r;
        
        % Update solution
        x = x + alpha * z;
        
        % Compute residual norm
        residual_norm = norm(r);
    end
end