N = 100;
relativeError = ones(2,100000);
relativeResidue = ones(2,100000);
count = 0;
iterations = ones(1,6);
k = 4;
[Q1, R ] = qr(randn(N));
[Q2, R ] = qr(randn(N));
D = diag(10 .^(k/N:k/N:k));
Aactual = Q1*D*Q1';
X = sort(Aactual')';
sort1 = X(:,1*floor(log2(N)));
sort2 = X(:,N-1*ceil(log2(N)));
sort1 = sort1*ones(1,N);
sort2 = sort2*ones(1,N);
f1 = 1.0*(Aactual<sort1) + 1.0*(Aactual>sort2);
A1 = f1.*Aactual ;


bactual = ones(N,1);
A = Aactual;
b = bactual;
x0 = zeros(N,1);
xnew = x0;
ynew = x0;
r = b - A*xnew;
r1 = b - A'*ynew;
p = r;
q = r1;
err=1;
it = 0;
eps = 0.00001;

while( err > eps)
    Ap = A*p;
    alpha = (r1'*r)/(q'*Ap);
    xtrue = xnew;
    xnew = xnew + alpha*p;
    ynew = ynew + alpha*q;
    rold = r;
    r1old = r1;
    r = r - alpha*Ap;
    r1 = r1 - alpha*A'*q;
    beta = (r1'*r)/(r1old'*rold);
    p = r + beta*p;
    q = r1 + beta*q;
    count = count + 1;
    err = norm((r)/norm(b));
    relativeError(1, count) = it;
    relativeError(2, count) = norm((xtrue - xnew)/norm(xnew));


    relativeResidue(1, count) = it;
    relativeResidue(2, count) = err;

    it = it+1;
end

disp(['No of iterations for non - preconditioned case = ' num2str(it)])
disp(['b - Ax error = ' num2str(norm(bactual-Aactual*xnew))])


setup.type = 'nofill';
setup.milu = 'row';
setup.droptol = 0.5;
[L,U] = ilu(sparse(A1),setup);


precond = L*U;
%precond = diag(1./diag(Aactual));
b = bactual;

x0 = zeros(N,1);
xnew = x0;
ynew = x0;
r = bactual - Aactual*xnew;
r1 = bactual - Aactual'*ynew;
p = precond*r;
q = precond*r1;
err=1;  
it = 0;
while( err > eps*10^-k)
    Ap = Aactual*p;
    alpha = (r1'*precond*r)/(q'*Ap);
    xtrue = xnew;
    xnew = xnew + alpha*p;
    ynew = ynew + alpha*q;
    rold = r;
    r1old = r1;
    r = r - alpha*Ap;
    r1 = r1 - alpha*Aactual'*q;
    beta = (r1'*precond*r)/(r1old'*precond*rold);
    p = precond*r + beta*p;
    q = precond*r1 + beta*q;
    count = count + 1;
    err = norm((r)/norm(b));
    relativeError(1, count) = it;
    relativeError(2, count) = norm((xtrue - xnew)/norm(xnew));


    relativeResidue(1, count) = it;
    relativeResidue(2, count) = err;

    it = it+1;
end
disp(['No of iterations for preconditioned case = ' num2str(it)])
disp(['b - Ax error = ' num2str(norm(bactual-Aactual*xnew))])