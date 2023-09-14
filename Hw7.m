%7.1

for n= [10, 20, 40]

    H=hilb(n);
    x=ones(n,1);
    b=H*x;

    alpha=1e-6;
    deltab=randn(n,1)*alpha;
    deltax=H\deltab;
   
    Kest= norm(deltax)/norm(x) * norm(b)/norm(deltab)
    cond(H,2)
   

end


 
n=[1:20];

for ii=n
    H=hilb(ii);
    xi=ones(ii,1);
    b=H*xi;
    x=H\b;
    norm(x-xi,2)/norm(x,2)
    ii
end



%7.4

clear all

%a

A=3*eye(6)-1*diag(ones(5,1),1)-1*diag(ones(5,1),-1);
A(1, 6) = -1;
A(2, 5) = -1;
A(5, 2) = -1;
A(6, 1) = -1;
A
b=[1 0 1 1 0 1]';

D=diag(diag(A));
E=-tril(A,-1);

Bgs=(D-E)\(D-E-A)
Bj=D\(D-A)

ggs=(D-E)\b
gj=D\b

x0=zeros(6,1);
tol=1e-8;
max_it= 1000;

max(abs(eig(Bgs)))<1
%A spd

max(abs(eig(Bj)))<1
%A 

[x, iter, incr] = stationary_method(Bgs, ggs, x0, tol, max_it);
iter
[x, iter, incr] = stationary_method(Bj, gj, x0, tol, max_it);
iter



