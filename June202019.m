%
n=10;
d=[1:n];
v=ones(n-1,1);
%1,0.5,1,05
for i=1:2:n-1
    v(i)=1;
end
for i=2:2:n-2
    v(i)=0.5;
end
A=diag(d)+diag(v,+1)+diag(v,-1);

D=diag(diag(A));
E=-tril(A,-1);
Bj=D\(D-A);
Bgs=(D-E)\(D-E-A);
max(abs(eig(Bj)))<1;
max(abs(eig(Bgs)))<1;

xex=ones(n,1);
b=A*b;

tol=1e-6;
maxit=500;

x0=zeros(n,1);

gj=D\b;
[x, iter, incr] =stationary_method(Bj,gj,x0,tol,maxit);
numel(iter)

ggs=(D-E)\b;
[x, iter, incr] =stationary_method(Bgs,ggs,x0,tol,maxit);
numel(iter)


P=diag(diag(A));
alpha=2/(max(abs(eig(inv(P)*A)))+min(abs(eig(inv(P)*A))));
[x, iter, incr] =prec_rich_method(A,b,P,alpha,x0,tol,maxit);
numel(iter)



%all three methods converge in just one iteration, meaning that the values
%of the preconditioned matrix P and alpha were properly chosen for the
%problem at hand: they yield good convergence properties




%
f1 =@(x) (x - 5).^2;
f2 =@(x) f1(x).^2;
a=0;
b=1;

x=linspace(a,b,1000);
y1=f1(x);
y2=f2(x);

figure;
plot(x,y1)
hold on;
plot(x,y2)
xlabel('x');
ylabel('f(x)');
grid on;

m=20;
composite_trapezoidal(f1,a,b,m)
61/3
composite_trapezoidal(f2,a,b,m)
2101/5
 

%both rules are quite accurate, in case of f1 the precision is higher since
%we are dealing with a polynomial of degree 2
%Since f1,f2 are at least C^2 we guarantee that both rules converge with an
%order of accuracy equal to 2

for i=1:5
    integer_m(i)=composite_trapezoidal(f2,a,b,2.^i)
end

err=abs(2101/5-integer_m);
p_M=-diff(log(err))/log(2);
p_M





