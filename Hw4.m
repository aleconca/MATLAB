A = [12 5 -8 -5; -4 -4 8 -6; 4 2 -3 0; 0 -1 2 -4]
b = [4 -6 3 -3]';

%a
[L, U]=lu(A);

y = forward_substitution(L,b)
x = backward_substitution(U,y)

A\b

%b
detA = prod(diag(L))*prod(diag(U))
det(A)

%c
c = [-22 -12 -1 -11]'

y = forward_substitution(L,c)
x = backward_substitution(U,y)

A\c
