%

n = 50;  % Size of the matrix
B=ones(n);
B(1,1:n)=3;
B(1:n,1)=3;
B(2:n, 2:n) = 5; 
for i=1:n
   B(i,i)=3;
end   

tol = 1e-6;
max_it = 500;
x =3*ones(50,1);
x0=zeros(50,1);

I=50*eye(50,50);
A=B+I;
b=A*x;

D=diag(diag(A));
E=-tril(A,-1);
Bj=D\(D-A);
gj=D\b;
max(abs(eig(Bj)))<1 %not convergent
[x, iter, incr] = stationary_method(Bj,gj,x0,tol,max_it);
iter

Bgs=(D-E)\(D-E-A);
ggs=(D-E)\b;
max(abs(eig(Bgs)))<1
[x, iter, incr] = stationary_method(Bgs,ggs,x0,tol,max_it);

P = diag(A); %Pj=D
alpha= 0,9154;
[x, iter, incr] =prec_rich_method(A,b,P,alpha,x0,tol,max_it);
iter



%
f =@(x) exp(-x+2).*log(x);
a=0.1;
b=4;

% Define the range of x values
x = linspace(a, b, 100);

% Evaluate the function
y = f(x);

% Plot the function
figure;
plot(x, y)
hold on

% Plot the zero of the function
zero_x = fzero(f, [a, b]);
zero_y = f(zero_x);
plot(zero_x, zero_y, 'ro')


% Add labels and title
xlabel('x')
ylabel('f(x)')
title('Plot of f(x) = e^{-x+2} * log(x) with its zero')

% Add gridlines
grid on



%
tol = 1e-6; 
Nmax = 250; 
x0 = 2.5;
df=@(x) -exp(-x+2).*log(x)+exp(-x+2)./(x);

x0=bisection(f,a,b,tol);
[x,x_iter]=newton(f,df,x0,tol,Nmax);
x
numel(x_iter)


%



