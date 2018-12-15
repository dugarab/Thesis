load('A_sparse_LU_sparse.mat');
% load('C:\Users\Abhishek\Desktop\IISc\Thesis\A_sparse_LU_sparse_GMRES.mat')
arr = [100 200 500];
K=9;
gain = zeros(K,3);
iterations_pre = zeros(K,3);
iterations = zeros(K,3);
counts = zeros(K,3);
j=1;
for k = 1:K
    j=1;
    for N = arr
        
        for i = 1:1000
            if(no_of_it_pre(k,i,j)==0)
                continue
            end
            gain(k,j) = gain(k,j) + 1.0*no_of_it(k,i,j)/no_of_it_pre(k,i,j) ;
            iterations_pre(k,j) = iterations_pre(k,j) + no_of_it_pre(k,i,j);
            iterations(k,j) = iterations(k,j) + no_of_it(k,i,j);
            counts(k,j) = counts(k,j) + 1.0;
           
        end
        gain(k,j) = gain(k,j)/counts(k,j);
        
        iterations(k,j) = iterations(k,j)/counts(k,j);
        iterations_pre(k,j) = iterations_pre(k,j)/counts(k,j);
                    
        j = j+1;
    end
end
figure(1)
subplot(3,1,1)
plot(1:K,gain(:,1)','-o')
title('Average gain for condition numbers (N=100)')
xlabel('condition number')
ylabel('Average gain')

subplot(3,1,2)
plot(1:K,gain(:,2)','-o')
title('Average gain for condition numbers (N=200)')
xlabel('condition number')
ylabel('Average gain')
subplot(3,1,3)
plot(1:K,gain(:,3)','-o')
title('Average gain for condition numbers (N=500)')
xlabel('condition number')
ylabel('Average gain')


figure(2)
subplot(3,1,1)
plot(1:K,iterations(:,1)','-o')
title('Average iterations for condition numbers (N=100)')
xlabel('condition number')
ylabel('Average iteration')
subplot(3,1,2)
plot(1:K,iterations(:,2)','-o')
title('Average iterations for condition numbers (N=200)')
xlabel('condition number')
ylabel('Average iteration')
subplot(3,1,3)
plot(1:K,iterations(:,3)','-o')
title('Average iterations for condition numbers (N=500)')
xlabel('condition number')
ylabel('Average iteration')

figure(3)
subplot(3,1,1)
plot(1:K,iterations_pre(:,1)','-o')
title('Average iterations for pre-conditioned system(N=100)')
xlabel('condition number')
ylabel('Average iteration')
subplot(3,1,2)
plot(1:K,iterations_pre(:,2)','-o')
title('Average iterations for pre-conditioned system(N=200)')
xlabel('condition number')
ylabel('Average iteration')
subplot(3,1,3)
plot(1:K,iterations_pre(:,3)','-o')
title('Average iterations for pre-conditioned system(N=500)')
xlabel('condition number')
ylabel('Average iteration')
        
        