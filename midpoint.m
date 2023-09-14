function int = midpoint(f, a, b)
	% MIDPOINT Computes the integral of a function f between a and b using 
	%   the midpoint quadrature rule.
	%
	%   Inputs  : f      = function handle of the function to be integrated
	%             a      = left limit of the interpolation interval
	%             b      = right limit of the interpolation interval
	%   Outputs : int = integral approximation

	int = (b-a)*f( (a+b)/2 );

end