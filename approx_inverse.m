N = 500;
rc = 1/10^4;
%A = sprand(N,N,N*log2(N)/(N*N),rc);
A = sprandsym(N,N*log(N)/(N*N),rc, 2);

M = sparse(eye(N));

%for i =1:N
%    if abs(A(i,i)) > 1e-4       
%        M(i,i) = 1/A(i,i);
%    end
%end
I = sparse(eye(N));
for i = 1:10
    R = I - A*M;
    G = A'*R;
    temp = A*G;
    alpha = (norm(G,'fro')/norm(temp,'fro'))^2;   
    M = M + alpha*G;

end

condest(A)
condest(A*M)
