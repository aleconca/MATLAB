clc
%9.1

% c)

f1 = @(x) x.^3;
f2 = @(x) x.^5;
a = 0;
b = 1;

% Midpoint rule: the number of subintervals is equal to the number of total nodes.

I1_m_1  = composite_midpoint(f1, a, b, 1)
I1_m_10 = composite_midpoint(f1, a, b, 10)

% The function to be integrated is a polynomial of order 3. Hence the rule,
% of degree 1, is not able to compute the integral exactly.
% The error decreases while increasing the number of nodes
% (i.e., increasing the number of subintervals)

% Trapezoidal rule: the number of subintervals is equal to the number of total nodes minus 1.

I1_t_2  = composite_trapezoidal(f1, a, b, 1)
I1_t_10 = composite_trapezoidal(f1, a, b, 9)

% Same behaviour as before, since the trapezoidal rule is of degree 1.

% d)

% Simpson rule: the number of subintervals is equal to (t - 1)/2,
% where t is the number of total nodes.

I1_s_3 = composite_simpson(f1, a, b, 1)
I1_s_7 = composite_simpson(f1, a, b, 3)

% The result does not depend on the number of nodes.
% The reason of this is the fact that Simpson rule has degree of exactness
% equal to 3 and the integral on any subinterval is computed exactly

I2_s_3 = composite_simpson(f2, a, b, 1)
I2_s_7 = composite_simpson(f2, a, b, 3)

% Instead, when considering the approximation of I_2,
% the integrand function is a polynomial of order higher than 3, hence the rule is not exact.

% e)

% Recall that the error for the composite trapezoidal rule is given in
% slide 6.

h_star = sqrt(12 * 1e-3 / ((b-a) * 20))
m_star = ceil((b-a) / h_star)

% ... and the number of nodes is m_star + 1.

% Check
I1_t_star = composite_trapezoidal(f1, a, b, m_star)
err_star = abs(1/4 - I1_t_star)

% Similarly for the composite Simpson rule:

h_star = ((180 * 1e-3) / ((b-a) * 120))^(1/4)
m_star = ceil((b-a) / (2*h_star))
% ... and the number of nodes is 2 m_star + 1.

err_star = abs(1/6 - I2_s_7)
pause

% f)

for (i = 0:9)
  m = 2^i;
  integr_m(i+1) = composite_midpoint(f2, a, b, m);
  integr_t(i+1) = composite_trapezoidal(f2, a, b, m);
  integr_s(i+1) = composite_simpson(f2, a, b, m);
end

err_m = abs(1/6 - integr_m);
err_t = abs(1/6 - integr_t);
err_s = abs(1/6 - integr_s);

p_m = -diff(log(err_m)) / log(2)%number of nodes x2 every time
p_t = -diff(log(err_t)) / log(2)
p_s = -diff(log(err_s)) / log(2)

% As expected, the error decreases quadratically in h for both midpoint
% and trapezoidal rules. For Simpson rule instead, convergence is of fourth order in h.






%9.2
f=@(x) (1-x.^2).^(1/2);
a=-1;
b=1;

for i=1:9
    comp_trap(i)=composite_trapezoidal(f,a,b,3.^i-1)
    comp_sim(i)=composite_simpson(f,a,b, (3.^i-1)./2)
end

errT=abs(pi/2-comp_trap)
errS=abs(pi/2-comp_sim)

p_T=- diff(log(errT)) /log(3)%number of nodes x3 every time
p_S=- diff(log(errS)) /log(3)






%9.3
% a)

weights_1 = 2/3*[2 -1 2];
nodes_1   = [-1/2 0 1/2];
f = @(x) x.^0; sum(weights_1 .* f(nodes_1))
f = @(x) x.^1; sum(weights_1 .* f(nodes_1))
f = @(x) x.^2; sum(weights_1 .* f(nodes_1))
f = @(x) x.^3; sum(weights_1 .* f(nodes_1))
f = @(x) x.^4; sum(weights_1 .* f(nodes_1)) % inexact
pause

weights_2 = 1/4*[1 3 3 1];
nodes_2   = [-1 -1/3 1/3 1];
f = @(x) x.^0; sum(weights_2 .* f(nodes_2))
f = @(x) x.^1; sum(weights_2 .* f(nodes_2))
f = @(x) x.^2; sum(weights_2 .* f(nodes_2))
f = @(x) x.^3; sum(weights_2 .* f(nodes_2))
f = @(x) x.^4; sum(weights_2 .* f(nodes_2)) % inexact

% Both formulae are exact up to degree 3.

% b)

% The integral can be rewritten as
% int_{-1}^{1} log(t+2) dt;
% we then apply the quadrature rule.

f = @(t) log(t+2)
Q1 = sum(weights_1 .* f(nodes_1))
Q2 = sum(weights_2 .* f(nodes_2))

F = @(x) x.*log(x) - x;
int_exact = F(3) - F(1)

err_Q1 = abs(int_exact - Q1)
err_Q2 = abs(int_exact - Q2)

% The second rule behaves better.






%9.4

% 1)
a = 0;
b = 1;
%alpha = 2;
alpha = 3/2;
%alpha = 5/2;
%alpha = 7/2;
%alpha = 2;
f = @(x) (x.^alpha);

% Use an increasing number of nodes in computing the integral
for (i = 0:9)
  m = 2^i;
  integr_m(i+1) = composite_midpoint(f, a, b, m);
  integr_t(i+1) = composite_trapezoidal(f, a, b, m);
  integr_s(i+1) = composite_simpson(f, a, b, m);
end

int_exact = quad(f, a, b, 1e-12);

err_m = abs(int_exact - integr_m);
err_t = abs(int_exact - integr_t);
err_s = abs(int_exact - integr_s);

p_m = -diff(log(err_m)) / log(2)
p_t = -diff(log(err_t)) / log(2)
p_s = -diff(log(err_s)) / log(2)

% Depending on the considered value of alpha the integrand function is characterized
% by a different level of regularity.
% With alpha = 1/2, all the methods show convergence of order 3/2 (= alpha + 1).
% The regularity of the function is too low for the methods to reach their
% maximum theoretical convergence order.
% With alpha = 3/2 the error decreases quadratically in h for both midpoint
% and trapezoidal rules. For the Simpson rule instead, convergence is of order 5/2 (= alpha + 1).
% Again regularity is not enough for the Simpson method to reach fourth order.
% Similar situation occurs with alpha = 5/2, with Simpson rule showing order 7/2 (= alpha + 1).
% Finally, with alpha = 7/2 all the methods show their theoretical orders.

