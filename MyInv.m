function [InvA] = MyInv(A)

    [L,U] = lu(A);
    %s = size(A);
    %n = s(1);
    n = length(A);
    
    for k = 1:n
        
        e = zeros(n,1);
        e(k) = 1;
        
        y = forward(L,e);
        InvA(:,k) = backward(U,y);
        
    end

end