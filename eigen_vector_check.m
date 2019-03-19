nums = [100,200,400,800];
cond_num = zeros(4,500);
vals = zeros(size(nums));
total_iters = 3;
j=0;
for N = nums
    j = j+1;
    for i = 1:total_iters
        A = randn(N);
        [K, ~] = eig(A);
        cond_num(j,i) = cond(K);
    end
    vals(j) = sum(cond_num(j,:))/total_iters;
end
plot(nums, vals);
xlabel('Size of matrix');
ylabel('Condition number of eigenVector matrix');
        
        