%
n=50;
A=eye(n)+(-0.2)*diag(ones(n-2,1),-2)+(-0.3)*diag(ones(n-2,1),+2)+0.1*diag(ones(n-3,1),-3)+0.1*diag(ones(n-3,1),+3);
for i=1:n-1
    det(A(1:i,1:i))~=0;%verified N.S. condition
end


[L,U,P]=lu(A)
P;%No pivoting because P is still the identity matrix->no rows were permutated

xex=[1:n];
b=A*b;

y=L\b;
x=U\y;

%
f=@(x) x.^2 - 3*x + 2;
a=0;
b=4;

x_plot=linspace(a,b,1000);
y_plot=f(x_plot);
x1=2;
x2=1;
y1=f(x1);
y2=f(x2);


figure;
plot(x_plot,y_plot);
xlabel('x');
ylabel('f(x)');
hold on;

plot(x1,y1, 'ro')
plot(x2,y2, 'ro')
grid on;


%
g1 =@(x) sqrt(3*x - 2); %Ostrowski theorem satisfied only foe zero x=2 of f
g2 =@(x) (x.^2 + 2)./3; %No fixed point

%maxit fixed point
tol=1e-6;
maxit=500;
x0=2;
[xi, x_iter] =fixed_point(g1,x0,tol,maxit);
numel(x_iter)
xi



%







