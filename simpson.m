function int = simpson(f, a, b)
	% SIMPSON Computes the integral of a function f between a and b using 
	%   the simpson quadrature rule.
	%
	%   Inputs  : f      = function handle of the function to be integrated
	%             a      = left limit of the interpolation interval
	%             b      = right limit of the interpolation interval
	%   Outputs : int = integral approximation

	int = (b-a)/6*(f(a) + 4*f((a+b)/2) + f(b));

end