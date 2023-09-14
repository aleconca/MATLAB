function [L,U,x] = thomas(A,b)
	%   [L,U] = lu_decomposition(A)
	%
	%   Inputs  : A = input tridiagonal matrix in sparse format
	%             b = input rhs vector
	%   Outputs : L = lower triangular matrix
	%             U = upper triangular matrix
	%             x = solution of the linear system

	n = length(b);

	c = spdiags(A, 1); % of the form c = [0, c_1, c_2, ... c_{n-1}]
	a = spdiags(A, 0); % of the form a = [a_1, a_2, ... a_n]
	e = spdiags(A,-1); % of the form e = [e_2, e_3, ... e_n, 0]

	y = zeros(n,1);
	x = zeros(n,1);
	alpha = zeros(n,1); % of the form alpha = [alpha_1, alpha_2, ... alpha_n]
	beta  = zeros(n,1); % of the form beta = [beta_2, beta_3, ... beta_n, 0]
	eyediag = ones(n,1);
		
	%% LU factorization %%
	alpha(1) = a(1);
	for i=2:n
	    beta(i-1) = e(i-1)/alpha(i-1);
	    alpha(i) = a(i) - beta(i-1)*c(i);
	end
	L = spdiags([beta eyediag], [-1 0], n, n);
	U = spdiags([alpha c], [0 1], n, n);
	
	%% Forward substitution %%
	y(1)=b(1);
	for i=2:n
	    y(i) = b(i) - beta(i-1)*y(i-1);
	end
	
	%% Backward substitution %%
	x(n) = y(n) / alpha(n);
	for i = n-1:-1:1
	    x(i) = (y(i) - c(i+1)*x(i+1))/alpha(i);
	end

end