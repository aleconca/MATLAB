%6.1
% We can invert the relation
% norm(deltax)/norm(x) <= norm(deltab)/norm(b)
% to get an estimate for Kp(A):
% Kp(A) >= norm(deltax)/norm(x) * norm(b)/norm(deltab)

for (n = [10 20 40])
	pert_mag = 1e-6; % Perturbation magnitude
	n_trials = 100;
	H = hilb(n);
	x = ones(n, 1);
	b = H*x;
	KH_est = 0;
	for (ii = 1:n_trials)
	  deltab = randn(n, 1) * pert_mag;
	  deltax = H \ deltab;
	  KH_est_ii = norm(deltax)/norm(x) * norm(b)/norm(deltab); 
	  KH_est = max(KH_est_ii, KH_est);
	end
	KH_est
end

cond( hilb(10) )
cond( hilb(20) )
cond( hilb(40) )

% The estimates obtained at the point a) are lower. This is expected, since
% the condition number represents the maximum 
% possible amplification of the relative error on data. With the first
% approach we might not have considered the worst condition, 
% hence the estimate is lower than the result with cond command.

n = [1:20];

for (ii = n)
  H = hilb(ii);
  x_ex = ones(ii, 1);
  b = H*x_ex;
  x = H\b;
  KH(ii) = cond(H);
  err(ii) = norm(x - x_ex) / norm(x);
end




%6.4
A = [1 2 -2;
     1 1  1;
     2 2  1];
b = A * [1 2 3]';

D  = diag(diag(A));
L  = tril(A, -1);
U  = triu(A, 1);
Bj = -D \ (L+U);
gj = D \ b;

x0 = zeros(3,1);

iter   = 0;
%tol    = 1e-5;
tol    = 1e-8;
maxit = 100;

[xj, iterj, incrj] = stationary_method(Bj, gj, x0, tol, maxit)

% Four iterations are performed in both cases, independently on the tolerance values we considered.
% Notice that the last increment is actually zero!

% The spectral radius is very small.

eig(Bj)
rhoBj  = max( abs( eig(Bj) ) )

% Indeed, computing by hand the eigenvalues we see that
% lambda_1 = lambda_2 = lambda_3 = 0, so the spectral radius is actually 0.

% Compute the cube of Bj
Bj^2
Bj^3        % At the third step we get the null matrix.

% Thus, for all k > 3
% Jacobi method has converged to the exact solution




%6.8,6.9
