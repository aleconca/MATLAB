function int = trapezoidal(f, a, b)
	% TRAPEZOIDAL Computes the integral of a function f between a and b using 
	%   the trapezoidal quadrature rule.
	%
	%   Inputs  : f      = function handle of the function to be integrated
	%             a      = left limit of the interpolation interval
	%             b      = right limit of the interpolation interval
	%   Outputs : int = integral approximation

	int = (b-a)/2*(f(a)+f(b));

end
