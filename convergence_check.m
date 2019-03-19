k = 3;
%t = -1:0.01:1;
%plot(t,chebyshevT(k,sym(t)));

mode = 'nofill';
N = 500;
avg = [];
iterations = [];
no_of_trials = 1;
val = (1+10^k)/(1-10^k);
w = val + (val^2 - 1)^0.5;
C = zeros(N,1);
for i = 1:N
    C(i) = 1/abs(0.5*(w^i + w^(-i)));
end

%C = 1./(abs(chebyshevT(1:N,vpa(val))));

for i = 1:no_of_trials
     offset = 0;

    [Q, ~] = qr(randn(N));
    temp = Q*(diag(10 .^( k/N:k/N:k)))*Q';
    A = temp;%*(diag(abs(randn(N,1))))*inv(temp);
    b = randn(N,1); b = b/norm(b);
    tol = 10^-3;
    %A = Aactual;
    m=N;
    [x1,FLAG,RELRES,iters1,e_res] = gmres(A,b,[],tol,N);
    %[x1, e_res, e_estimator, iters1, x_ans] = gmres_custom(A,b,zeros(N,1),m,tol);
    r_0 = norm(b - A*x0);
    alpha = 10^(k/N);
    beta = 10^(k);
    %iters1 = iters1-1;
    iters1 = iters1(2);
    
    plot(1:iters1, e_res(1:iters1), 1:iters1, C(1:iters1)) 
    legend('Actual', 'Calculated');
    %plot(1:no_of_trials, avg);
    title('Ratio of calculated vs actual rel residue (cond = 10^3, N = 500)')
    ylabel('Ratio');
    xlabel('Iteration number');

    
    %tic();

    
end
%plot(1:iters1, e_res, 1:iters1, C(1:iters1)) 
%legend('calculated', 'actual');
%plot(1:no_of_trials, avg);
%title('Ratio of calculated vs actual rel error (cond = 10^4, N = 500)')
%ylabel('Ratio');
%xlabel('Iteration number');
%Avg_val = sum(avg)/no_of_trials
