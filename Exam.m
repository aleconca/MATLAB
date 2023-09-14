%6-02-23
A=[ 2, 0, 2; 
    0, 12, 0; 
    2, 0, 26];

D = diag(diag(A))
E=-tril(A,-1);
Bj=D\(D-A);
Bgs=(D-E)\(D-E-A);

rhoj=max(abs(eig(Bj))) <1; % sufficient conditions not satisfied
rhogs=max(abs(eig(Bgs)))<1;


A1=[12, -3, 10;
    1, 6, 0; 
    8, 1, -11];

D = diag(diag(A1));
E=-tril(A1,-1);
Bj1=D\(D-A1);

Bgs1=(D-E)\(D-E-A1);

rhoj1=max(abs(eig(Bj1)))<1;
rhogs1=max(abs(eig(Bgs1)))<1;


%prec richardson scheme
n=15;
A = 2*diag(ones(n,1))-1*diag(ones(n-1,1),-1)-1*diag(ones(n-1,1),1)
b=A*ones(n,1);
x0=zeros(n,1);
tol=1e-5;
maxit=1000;
P=diag(diag(A))%475 ierations
P1=diag(diag(A))+1*diag(ones(n-1,1),1)%not symmetric

alpha=2/(max(eig(inv(P)*A))+min(eig(inv(P)*A)));

[x, iter, incr] = prec_rich_method(A, b, P1, alpha, x0, tol, maxit)



%fixed point
phi1=@(x) 0.5 * (x-x^3);
phi2=@(x) 0.5 * log(1-x);
alpha = 0;
tol =1e-6;
x0=0.5;
maxit=200;

[xi1, x_iter1]=fixed_point(phi1, x0, tol, maxit)%faster
numel(x_iter1)
[xi2, x_iter2]=fixed_point(phi2, x0, tol, maxit)
numel(x_iter2)


%

f=@(x) abs(x.^1/3);
a=-2;
b=2;

%uniform intervals

%piecewise linear interpolant
n_vect=[2,4,8,16,32];

for i = 1:numel(n_vect)
 x_nodes=linspace(a,b,n_vect(i)+1);
 y_nodes=f(x_nodes);

 x_plot=linspace(a,b,1000);
 p_plot=interp1(x_nodes,y_nodes,x_plot);

 figure;
 plot(x_plot,f(x_plot),x_plot,p_plot, x_nodes,y_nodes)%plot multiple data
 err(i)=max(abs(f(x_plot)-y_plot))
end

figure
loglog(H, err, 'bx-', 'LineWidth', 2, 'MarkerSize', 8)
hold on, box on
loglog(H, H.^(2), 'g-', 'LineWidth', 2)
axis([1e-1 1e1 1e-3 1e2])
set(gca,'FontSize',16)
set(gca,'LineWidth',1.5)
xlabel('h','FontSize',16)
ylabel('error','FontSize',16)
legend('error', 'O(h^2)', 'Location', 'ne');

p=-diff(log(err)) / log(2) %x2 each time


for i=0:8
x_i_hat = -cos(pi*i./n);
x_i=(a+b)./2 + (b-a)*x_i_hat;

%x_nodes=linspace(a,b,8+1);
y_nodes=f(x_i);
p=polyfit(x_i,y_nodes,8);

x_plot=linspace(a,b,1000);
%y_plot=f(x_plot);
p_plot=polyval(p,x_plot);

figure
plot(x_plot,f(x_plot),x_plot,p_plot, x_i,y_nodes)%plot multiple data

end


%




