%2.1 Bisection
epsilon = 0.01;

tol = 1e-6;
h = @(x) cot(x) - (x.^2 - 1)./(2*x);

rootfinding_function_plot(h, 0+epsilon, pi-epsilon, true);
rootfinding_function_plot(h, pi+epsilon, 2*pi-epsilon, false);
rootfinding_function_plot(h, 2*pi+epsilon, 3*pi-epsilon, false);

I = [6.5 7] % cannot use [6, 7] since h is discontinuous there
[xi, x] = bisection(h, I(1), I(2), tol);
xi
iter = numel(x)%returns the number of elements in x

% Since the exact value of the roots is NOT known, we perform
% again the bisection method with a very small tolerance, and consider
% the final approximation as the exact root of the equation.

figure

% First root
[xiex, xex] = bisection(h, I(1), I(2), 1e-10*tol);
xiex
err = abs(x - xiex);

%subplot(1,2,1)%grid dimension, se aumenti il terzo valore schiacci il grafico lungo x
semilogy(err, 'LineWidth',2)
%box on, hold on
%semilogy(tol*ones(iter,1), 'r--', 'LineWidth',2)
xlim([0 iter+1])%x axis range
%set(gca,'FontSize',16)
xlabel('Iteration','FontSize',16)
ylabel('Error','FontSize',16)


%2.2
f = @(x) exp(-(x-2).^2) + exp(x-4) - 1;
rootfinding_function_plot(f, 0, 5);

% There are three solutions in the interval [1, 5].

tol = 1e-6;

I1 = [0 5];
[xi1, x_iter1] = bisection(f, I1(1), I1(2), tol);
xi1

I2 = [1 6];
[xi2, x_iter2] = bisection(f, I2(1), I2(2), tol);
xi2

I3 = [1 5];
[xi3, x_iter3] = bisection(f, I3(1), I3(2), tol);
xi3

I4 = [-1 5];
[xi4, x_iter4] = bisection(f, I4(1), I4(2), tol);
xi4


% In the first and fourth interval the bisection method converges to the root with
% minimum magnitude; similarly, in the second and third interval the bisection method
% converges to the root with maximum magnitude. The second root cannot be reached.
% 
% This is caused by the iterative procedure that halves the interval: at the 
% very first iteration the half-interval that contains the second root is discarded.

bisection_plot(f, I2(1), I2(2), tol);
bisection_plot(f, I3(1), I3(2), tol);

% This example shows that the choice of the initial interval is important.


%2.3 Newton

f = @(x) tan(x) - 2*x;
epsilon = 0.01;
rootfinding_function_plot(f, 0+epsilon, pi/2-epsilon, true);

% From the plot we choose the interval I = [1, 1.5] that ensures the validity of the Bolzano's
% theorem, and we can use e.g. the bisection method to find x*.

xstar = bisection(f, 1, 1.5, 1e-15);


% The Newton method for the approximation of $\sin(x) = 0$ is 
% 
% x{(k+1)} = x{(k)} - {sin(x{(k)})}/{cos(x{(k)})} = x{(k)} - tan(x{(k)}).
%
% Therefore, if x{(0)} = x*
% 
% x{(1)} = x{(0)} - tan(x{(0)}) =  x* - tan(x*) =  x* - 2x* = -x* = -x{(0)}\\
% x{(2)} = x{(1)} - tan(x{(1)}) = -x* + tan(x*) = -x* + 2x* =  x* =  x{(0)}
% 
% The method is describing the orbit {x*, - x*}, and thus cannot converge!


g = @(x) sin(x);
dg = @(x) cos(x);

[xi, x_iter] = newton25(g, dg, xstar, 1e-6, 25);
x_iter

% The first 10 iterations of the Newton method seem to describe the orbit {x*, - x*}, but
% after the first few iterations some numerical error arises and this makes the algorithm move from the orbit and then to converge to the solution.


