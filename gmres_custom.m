function [x_final, e_residue, e_estimator, iters, x_ans] = gmres_custom( A, b, x_start, max_iterations, threshold, L, U)
  n = length(A);
  m = max_iterations;
  if ~exist('L', 'var')
      L = eye(n);
      U = eye(n);
  end
  x = zeros(n);
  x_ans = A\b;
  %use x as the initial vector
  r=L\(b-(A*(U\x_start)));

  b_norm = norm(b);
  error = norm(r)/b_norm;

  %initialize the 1D vectors
  sn = zeros(m,1);
  cs = zeros(m,1);
  e1 = zeros(n,1);
  e1(1) = 1;
  e=[error];
  r_norm=norm(r);
  Q(:,1) = r/r_norm;
  beta = r_norm*e1;
  iters = 0;
  for k = 1:m              
  
    if ~exist('L', 'var')
        [H(1:k+1,k) Q(:,k+1)] = arnoldi(A, Q, k);
    %run arnoldi
    else
        [H(1:k+1,k) Q(:,k+1)] = arnoldi(A, Q, k, L, U);
    end
    
    %eliminate the last element in H ith row and update the rotation matrix
    [H(1:k+1,k) cs(k) sn(k)] = apply_givens_rotation(H(1:k+1,k), cs, sn, k);
    
    %update the residual vector
    beta(k+1) = -sn(k)*beta(k);
    beta(k)   = cs(k)*beta(k);
    error  = abs(beta(k+1)) / b_norm;
    y = H(1:k,1:k) \ beta(1:k);
    x(:,iters+1) = x_start + U\(Q(:,1:k)*y);
    
    %save the error
    e=[e; error];
    iters = iters + 1;
    if ( error <= threshold)
      break;
    end
  end
  
  e_residue = zeros(1,iters-1);
  e_estimator = zeros(1,iters-1);
  norm_x = norm(x_ans);
  norm_b = norm(b);
  for iter = 1:iters-1
      e_residue(iter) = norm(b - A*x(:,iter))/norm_b;
      e_estimator(iter) = norm(x_ans-x(:,iter))/norm_x;
      
  end
  
  x_final = x(:,iters);
  %calculate the result
   
  
end




