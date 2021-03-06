simple_testing = 1;
method_name = 'gmres';

% %For GMRES sparse 
% condition_number_range = 6:10;
% runs_per_matrix_size = 1:1000;
% matrix_sizes = [100 200 500];

%For simple testing
simple_testing = 1;
condition_number_range = 1:8;
runs_per_matrix_size = 1:10;
matrix_sizes = [500];

normal = true;
mode = 'nofill';
no_of_it_pre = zeros(9,10000,3);
no_of_it = zeros(9,10000,3);

for k = condition_number_range
    for i = runs_per_matrix_size
        j =1;
        for N = matrix_sizes

            count = 0;
            iterations = ones(1,6);

            [Q1, R ] = qr(randn(N));
            [Q2, R ] = qr(randn(N));
            D = diag(10 .^(k/N:k/N:k));
            if normal
                Aactual = Q1*D*Q1';

            else
                Aactual = Q1*D*Q2;
            end
            
            if simple_testing == 1
                %when doing crap do it here
                %Aactual = Q1*diag(abs(randn(N,1)))*Q1';
                %Aactual = Q1*D*Q2';
                %Aactual = randn(N);
                %Aactual = sprand(N,N,(4*N*log2(N))/(N*N) ) + diag(randn(N,1));
            end

            X = sort(Aactual')';
            sort1 = X(:,2*floor(log2(N)));
            sort2 = X(:,N-2*ceil(log2(N)));
            sort1 = sort1*ones(1,N);
            sort2 = sort2*ones(1,N);
            f1 = 1.0*(Aactual<sort1) + 1.0*(Aactual>sort2) ;
            A1 = f1.*Aactual + diag(diag(Aactual));
            if simple_testing == 1
                %A1 = sprand(N,N,(4*N*log2(N))/(N*N) ) + diag(randn(N,1));
            end
            [L,U] = ilu(sparse(A1),struct('type',mode));
%            temp = A1;
%            A1 = Aactual;
            k1 = round(log10(cond(A1)));
%             k1=k;
            if simple_testing == 1
                
                disp(['Norm of (A - A_T)/norma(A) = ' num2str(norm(A1-A1')/norm(A1))])
                disp(['conition number of A = ' num2str(cond(A1))])
            end
            
            % setup.type = 'nofill';
            % setup.milu = 'row';
            % setup.droptol = 0.5;
            % [L,U] = ilu(sparse(A1),setup);
            if k1 < 1 && simple_testing == 0
                continue
            end

            bactual = randn(N,1);
            A = A1;
            b = bactual;

            tol = 0.0001*(10^-k1);
            maxit = N;
            m = N;
            if strcmp(method_name, 'gmres')
                [x0,fl0,rr0,it0,rv0] = gmres(A,b,m,tol,maxit);
                no_of_it(k1,i,j) = (it0(1)-1)*m+it0(2);
            elseif strcmp(method_name, 'bicg')
                [x0,fl0,rr0,it0,rv0] = bicg(A,b,tol,maxit);
                no_of_it(k1,i,j) = it0;
            end
            
            if simple_testing == 1
                disp(['No of iterations for non - preconditioned case = ' num2str(it0)])
                disp(['b - Ax error = ' num2str(norm(b-A*x0))])
            end
            
            if strcmp(method_name, 'gmres')
                
                [x1, e_res, e_estimator, iters1] = gmres_custom(A,b,zeros(N,1), maxit, tol,L,U);
                no_of_it_pre(k1,i,j) = iters1;
            elseif strcmp(method_name, 'bicg')
                [x1,fl1,rr1,it1,rv1] = bicg(A,b,tol,maxit,L,U);
                no_of_it_pre(k1,i,j) = it1;
            end
           
            
            if simple_testing == 1
                disp(['No of iterations for preconditioned case = ' num2str(it1)])
                disp(['b - Ax error = ' num2str(norm(b-A*x1))])
            end
            
            if i==1
                disp(['k = ' num2str(k) ' complete'])
            end
            
            j= j+1;
            
            %rank(full(L*U))
            
        end
    end
end
