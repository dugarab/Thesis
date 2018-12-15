N = 1000;
mode = 'nofill';
k=2.5;


[Q, ~] = qr(randn(N));
[Q2, ~] = qr(randn(N));
temp = Q*(diag(10 .^(k/N:k/N:k)))*Q2;
Aactual = temp*(diag(abs(randn(N,1))))*inv(temp);



% normal = true;
% [Q1, ~ ] = qr(randn(N));
% [Q2, ~ ] = qr(randn(N));
% D = diag(10 .^(k/N:k/N:k));
% if normal
%     Aactual = Q1*D*Q1';
% 
% else
%     Aactual = Q1*D*Q2;
% end

temp = abs(Aactual);
b = randn(N,1);
X = sort(temp')';
sort1 = X(:,2*floor(log2(N)));
sort2 = X(:,N-2*ceil(log2(N)));
sort1 = sort1*ones(1,N);
sort2 = sort2*ones(1,N);
f1 = 1.0*(temp<sort1) + 1.0*(temp>sort2) ;
A1 = (f1.*Aactual + diag(diag(Aactual)));
[L,U] = ilu(sparse(A1),struct('type',mode));


%A = sparse(A1);
A = Aactual;
tic();
[X,FLAG,RELRES,ITER] = gmres(A,b,[],10^-5,N);
toc();
disp(['No of iterations for non - preconditioned case = ' num2str(ITER)])
disp(['b - Ax error = ' num2str(norm(b-A*X)/norm(b))])

tic();
x = A\b;
toc();

tic();
[x1, e1, iters1] = gmres_custom(A,b, zeros(N,1),N, 10^-5,L,U);
toc();
disp(['No of iterations for preconditioned case = ' num2str(iters1)])
disp(['b - Ax error = ' num2str(norm(b-A*x1)/norm(b))])