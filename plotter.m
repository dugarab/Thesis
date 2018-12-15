load('A_dense_LU_sparse')
%load('C:\Users\Abhishek\Desktop\IISc\Thesis\A_dense_sparse_LU_GMRES.mat');
arr = [100 200 500];
gain = zeros(5,3);
iterations_pre = zeros(5,3);
iterations = zeros(5,3);
j=1;
for k = 2:6
    j=1;
    for N = arr
        
        for i = 1:1000
            
            gain(k-1,j) = gain(k-1,j) + 1.0*no_of_it_pre(k-1,i,j)/no_of_it(k-1,i,j) ;
            iterations_pre(k-1,j) = iterations_pre(k-1,j) + no_of_it(k-1,i,j)/1000.0;
            iterations(k-1,j) = iterations(k-1,j) + no_of_it_pre(k-1,i,j)/1000.0;
        end
        gain(k-1,j) = gain(k-1,j)/1000.0;
        %iterations(k-1,j) = iterations(k-1,j)/1000.0;
        %iterations_pre(k-1,j) = iterations_pre(k-1,j)/1000.0;
                    
        j = j+1;
    end
end
figure(1)
subplot(3,1,1)
plot(2:6,gain(:,1)','-o')
title('Average gain for condition numbers (N=100)')
xlabel('condition number')
ylabel('Average gain')

subplot(3,1,2)
plot(2:6,gain(:,2)','-o')
title('Average gain for condition numbers (N=200)')
xlabel('condition number')
ylabel('Average gain')
subplot(3,1,3)
plot(2:6,gain(:,3)','-o')
title('Average gain for condition numbers (N=500)')
xlabel('condition number')
ylabel('Average gain')


figure(2)
subplot(3,1,1)
plot(2:6,iterations(:,1)','-o')
title('Average iterations for condition numbers (N=100)')
xlabel('condition number')
ylabel('Average iteration')
subplot(3,1,2)
plot(2:6,iterations(:,2)','-o')
title('Average iterations for condition numbers (N=200)')
xlabel('condition number')
ylabel('Average iteration')
subplot(3,1,3)
plot(2:6,iterations(:,3)','-o')
title('Average iterations for condition numbers (N=500)')
xlabel('condition number')
ylabel('Average iteration')

figure(3)
subplot(3,1,1)
plot(2:6,iterations_pre(:,1)','-o')
title('Average iterations for pre-conditioned system(N=100)')
xlabel('condition number')
ylabel('Average iteration')
subplot(3,1,2)
plot(2:6,iterations_pre(:,2)','-o')
title('Average iterations for pre-conditioned system(N=200)')
xlabel('condition number')
ylabel('Average iteration')
subplot(3,1,3)
plot(2:6,iterations_pre(:,3)','-o')
title('Average iterations for pre-conditioned system(N=500)')
xlabel('condition number')
ylabel('Average iteration')
        
        