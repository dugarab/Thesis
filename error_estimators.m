
simple_testing = 0;
mode = 'nofill';
condition_number_range = [11];
runs_per_matrix_size = 1:500;
matrix_sizes = [200,500,1000];
no_of_it_pre_est = zeros(10000,3);
no_of_it_pre = zeros(10000,3);

for k = condition_number_range
    for i = runs_per_matrix_size
        j =1;
        for N = matrix_sizes

            [Q, ~] = qr(randn(N));
            temp = Q*(diag(10 .^(k/N:k/N:k)))*Q';
            Aactual = temp;%*(diag(abs(randn(N,1))))*inv(temp);
            b = randn(N,1);
            A1 = sparsiofy(Aactual,'mid');
            [L,U] = ilu(sparse(A1),struct('type',mode));

            A = sparse(A1);
            tol = 10^-8;
            %A = Aactual;

            %tic();
            [x1, e_res, e_estimator, iters1] = gmres_custom(A,b, zeros(N,1),N, tol,L,U);
            %toc();
            if simple_testing ==1 
                disp(['No of iterations for preconditioned case = ' num2str(iters1)])
                disp(['b - Ax error = ' num2str(norm(b-A*x1)/norm(b))])
            end
            iters2 = iters1;
            for it = 1:iters1-1
                if e_estimator(it) < tol^0.5
                    iters2 = it;
                    break;
                end
            end
            for it = 1:iters1-1
                if e_res(it) < tol^0.5
                    iters1 = it;
                    break;
                end
            end
            no_of_it_pre(i,j) = iters1;
            no_of_it_pre_est(i,j) = iters2;
            j=j+1;
            
        end
        if mod(i,10) == 0
            disp(['Run number: ' num2str(i) ' done!']);
        end
    end
end
semilogy(1:iters1-1, e_res, 1:iters1-1, e_estimator)
legend('residue', 'estimator')