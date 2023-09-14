%

f=@(x) exp(-abs(x).^2);
a=-4;
b=4;

figure;
x=linspace(a,b,1000);
y=f(x);
plot(x,y)
xlabel('x')
ylabel('f(x)')
title('Plot of f(x)')

n=9;%nodes
x_nodes=linspace(a,b,n);
y_nodes=f(x_nodes);
p=polyfit(x_nodes,y_nodes,n);

x_plot=linspace(a,b,1000);
%y_plot
interp=polyval(p,x_plot);

plot(x_nodes, y_nodes, x_plot,interp)

norm(f(x_plot)-interp,'inf')




n=100;
d=2;
x_nodes=linspace(a,b,n);
y_nodes=f(x_nodes);
p=polyfit(x_nodes,y_nodes,d);

x_plot=linspace(a,b,1000);
leastsquares=polyval(p,x_plot);

plot(x_plot,leastsquares)


norm(f(x_plot)-leastsquares,'inf')


%By comparing the Linfinity-norm values of the interpolation error and the least 
% square interpolation error, you can assess the accuracy 
% and effectiveness of the two methods for approximating the function f(x) in the given interval Ω.



%


