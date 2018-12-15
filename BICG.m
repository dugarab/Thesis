clear all;
N=100;
relativeError = ones(2,100000);
relativeResidue = ones(2,100000);

iterations = ones(1,6);
count = 0;

for k = 1:5
    [Q1, R ] = qr(randn(N));
    [Q2, R ] = qr(randn(N));
    D = diag(10 .^(k/N:k/N:k));
    A = Q1*D*Q2;
    b = ones(N,1);
    xtrue = A ^ (-1) *b;
    %M = [0 A ; A' 0];
    %c = [b ; b];

    x0 = ones(N,1);
    xnew = x0;
    ynew = x0;
    r = b - A*xnew;
    r1 = b - A'*ynew;
    p = r;
    q = r1;

    it = 0;
    while( norm(r) > 0.0001)
        alpha = (r1'*r)/(q'*A*p);
        xnew = xnew + alpha*p;
        ynew = ynew + alpha*q;
        rold = r;
        r1old = r1;
        r = r - alpha*A*p;
        r1 = r1 - alpha*A'*q;
        beta = (r1'*r)/(r1old'*rold);
        p = r + beta*p;
        q = r1 + beta*q;
        count = count + 1;
        relativeError(1, count) = it;
        relativeError(2, count) = norm((xtrue - xnew)/norm(xtrue));
        

        relativeResidue(1, count) = it;
        relativeResidue(2, count) = norm((r)/norm(b));
        
        it = it+1;
    end
    iterations(k) = count;
end

figure(1)
plot(relativeError(1,1:iterations(1)), log(relativeError(2,1:iterations(1))), ...
    relativeError(1,iterations(1)+1:iterations(2)), log(relativeError(2,iterations(1)+1:iterations(2))),...
    relativeError(1,iterations(2)+1:iterations(3)), log(relativeError(2,iterations(2)+1:iterations(3))),...
    relativeError(1,iterations(3)+1:iterations(4)), log(relativeError(2,iterations(3)+1:iterations(4))),...
    relativeError(1,iterations(4)+1:iterations(5)), log(relativeError(2,iterations(4)+1:iterations(5))));
legend('k=1','k=2','k=3','k=4','k=5','location' , 'northeast');
xlabel('No of Iterations')
ylabel('Relative Error')
title('Bi-CG Relative Error (log scale)');

figure(2)
plot(relativeResidue(1,1:iterations(1)), log(relativeResidue(2,1:iterations(1))), ...
    relativeResidue(1,iterations(1)+1:iterations(2)), log(relativeResidue(2,iterations(1)+1:iterations(2))),...
    relativeResidue(1,iterations(2)+1:iterations(3)), log(relativeResidue(2,iterations(2)+1:iterations(3))),...
    relativeResidue(1,iterations(3)+1:iterations(4)), log(relativeResidue(2,iterations(3)+1:iterations(4))),...
    relativeResidue(1,iterations(4)+1:iterations(5)), log(relativeResidue(2,iterations(4)+1:iterations(5))));
legend('k=1','k=2','k=3','k=4','k=5','location' , 'northeast');
xlabel('No of Iterations')
ylabel('Relative Residue')
title('Bi-CG Relative Residue (log scale)');
figure(2)
plot(1:it, relativeResidue);