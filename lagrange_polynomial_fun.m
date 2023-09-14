% NA LAB 2018-2019
% Authors A.G.Mauri
% Lagrange Polynomial function

function y = lagrange_polynomial_fun(xi,i,x)

n = size(xi,2);

z = 1 ; % because is a products

for j=1:n
    if (i~=j)
    z=z.*(x-xi(j))/(xi(i)-xi(j));
    end
end
   
y=z ;

end
