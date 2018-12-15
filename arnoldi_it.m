function [H,Q] = arnoldi_it(A,b,N)
    [~,row] = size(A);
    col = N;
    H = zeros(col+1,col);
    Q = zeros(row,col+1);
    q1 = b/norm(b);
    
    for n = 1:N
       Q(:,n) = q1;
       v = A*q1;
       
       for j = 1:n
          H(j,n) = Q(:,j)'*v;
          v = v -  H(j,n)*Q(:,j);          
       end
       
       H(n+1,n) = norm(v);
       q1 = v/H(n+1,n);       
       
    end
    Q(:,N+1) = q1;
end