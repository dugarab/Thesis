for cond_number = [3]
    for tol = [10^-4, 10^-5, 10^-6]


        testing = 0;
        N = 1000;
        rc = 1/10^cond_number;
        matrix_type = 2;
        total_iters = 1000;
        b = randn(N,1);
        b = b/norm(b);
        residue_cg1 = zeros(1,total_iters);
        residue_cg2 = zeros(1,total_iters);% zeros(200,1);
        total_iters_cg = zeros(1,total_iters);
        total_iters_gmres = zeros(1,total_iters);
        condition_num_gain = zeros(1,total_iters);%zeros(200,1);
        iteration_gain_cg = zeros(1,total_iters);%zeros(200,1);
        iteration_gain_gmres = zeros(1,total_iters);%zeros(200,1);

        %cond_gain = zeros(
        for i = 1:total_iters
            no_error = 1;
            while no_error == 1
                try        
                    R = sprandsym(N,N*log(N)/(N*N),rc, matrix_type);
                    L = ichol(R);
                    no_error=0;
                catch
                    if testing==1
                        disp('Caught Error');
                    end
                    no_error = 1;

                end
            end
            x_actual = R\b;
            [x1,flag1,relres1,iter1] = pcg(R,b,tol,N,[]);
            [x2,flag2,relres2,iter2] = pcg(R,b,tol,N,L,L');
            L_inv = inv(L);
            condition_gain = condest(L_inv*R*L_inv')/condest(R);
            condition_num_gain(i) = condition_gain;
            if testing ==1
                disp(['Condition gain reduction: ' num2str(condition_gain)]);
                disp(['Iteration gain cg: ' num2str(iter1/iter2)]);
                disp(['Residue: ' num2str((x2-x_actual)'*R*(x2-x_actual)/(x_actual'*R*x_actual))]);
                disp([iter1 iter2]);
            end
            residue_cg1(i) = (x2-x_actual)'*R*(x2-x_actual)/(x_actual'*R*x_actual);
            residue_cg2(i) = relres2;
            iteration_gain_cg(i) = iter1/iter2;
            total_iters_cg(i) = iter1;


            [x_final1, e_residue1, e_estimator1, iters1, x_ans1] = gmres_custom(R,b,zeros(N,1),N,tol);
            [x_final2, e_residue2, e_estimator2, iters2, x_ans2] = gmres_custom(R,b,zeros(N,1),N,tol,L,L');
            iteration_gain_gmres(i) = iters1/iters2;
            total_iters_gmres(i) = iters1;

            if testing ==1    
                disp(['Iteration gain gmres: ' num2str(iters1/iters2)]);
                disp(['Residue: ' num2str(norm((b-R*x_final1)))]);% ' ' num2str(relres2_g)]);
                disp([iters1 iters2]);
            end
        end
        filename = ['cg_gmres_cholesky_cond' num2str(cond_number) '_N' num2str(N) '_tol' num2str(-log10(tol))];
        save(filename);
    end
end
%cond_gain = [

