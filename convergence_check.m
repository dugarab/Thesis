k = 4;
%t = -1:0.01:1;
%plot(t,chebyshevT(k,sym(t)));

mode = 'nofill';
sizes = [100,200,500,1000];
iterations = [];
no_of_trials = 1000;
val = (1+10^k)/(1-10^k);
w = val + (val^2 - 1)^0.5;
iteration_numbers = zeros(length(sizes), no_of_trials);

%C = 1./(abs(chebyshevT(1:N,vpa(val))));
figure(1)
size_number = 1;
for N = sizes
    C = zeros(N,1);

    avg = zeros(1,no_of_trials);
    for i = 1:N
        C(i) = 1/abs(0.5*(w^i + w^(-i)));
    end

    for i = 1:no_of_trials
         offset = 0;

        %[Q, ~] = qr(randn(N));
        %temp = Q*(diag(10 .^( k/N:k/N:k)))*Q';
        %A = temp;%*(diag(abs(randn(N,1))))*inv(temp);
        A = sprandsym(N,N*log(N)/(N*N),1/(10^k), 2);
        b = randn(N,1); b = b/norm(b);
        tol = 10^-3;
        %A = Aactual;
        m=N;
        [x1,iters1,e_res] = cg_custom(A,b,tol);
        %[x1, e_res, e_estimator, iters1, x_ans] = gmres_custom(A,b,zeros(N,1),m,tol);
        %r_0 = norm(b - A*x0);
        alpha = 10^(k/N);
        beta = 10^(k);
        %iters1 = iters1-1;
        %iters1 = iters1(2);
        iteration_numbers(size_number, i) = iters1;

        %plot(1:iters1, e_res(1:iters1), 1:iters1, C(1:iters1)) 
        %legend('Actual', 'Calculated');
        %plot(1:no_of_trials, avg);
        %title('Ratio of calculated vs actual rel residue (cond = 10^3, N = 500)')
        %ylabel('Ratio');
        %xlabel('Iteration number');
        avg(i) = sum(abs(C(1:iters1)./e_res(1:iters1)'))/iters1;

        %tic();


    end
    disp(['N = ' num2str(N)]);
    subplot(length(sizes),1,size_number);
    plot(1:no_of_trials, avg);
    %legend('calculated', 'actual');
    %plot(1:no_of_trials, avg);
    %title('Difference between actual and calculated residue (cond = 10^' str2num(k) ', N = ' str2num(N) ')')
    %ylabel('Ratio');
    %xlabel('Iteration number');
    xlabel('Sample number');
    ylabel('Average Ratio');
    title(['Theoretical/Observed residue (cond = 10^' num2str(k) ', N = ' num2str(N) ')']);
    Avg_val = sum(avg)/no_of_trials
    size_number = size_number + 1;

end

figure(2)
size_number = 1;
for N = sizes
    
    subplot(length(sizes),1,size_number);
    plot(1:no_of_trials, iteration_numbers(size_number,:));
    temp = sum(iteration_numbers(size_number,:))/no_of_trials;
    xlabel('Sample number');
    ylabel('Average Ratio');
    title(['Iteration count (cond = 10^' num2str(k) ', N = ' num2str(N) ') avg = ' num2str(temp)]);
    Avg_val = sum(avg)/no_of_trials
    size_number = size_number + 1;

end
