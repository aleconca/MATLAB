function int = composite_trapezoidal(f, a, b, m)
	% COMPOSITE_TRAPEZOIDAL Computes the integral of a function f between a and b using 
	%   composite trapezoidal quadrature rule on m subintervals.
	%
	%   Inputs  : f      = function handle of the function to be integrated
	%             a      = left limit of the interpolation interval
	%             b      = right limit of the interpolation interval
	%             m      = number of subintervals
	%   Outputs : int = integral approximation

	x = linspace(a, b, m+1);
	int = 0;

	for (i = 1:m)
		int = int + trapezoidal(f, x(i), x(i+1));
	end

end
