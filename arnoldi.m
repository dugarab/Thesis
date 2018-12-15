
%----------------------------------------------------%
%                  Arnoldi Function                  %
%----------------------------------------------------%
function [h, q] = arnoldi(A, Q, k, L, U)
  if ~exist('L', 'var')
      q = A*Q(:,k);
  else
      q = L\(A*(U\Q(:,k)));
  end
  
  for i = 1:k
    h(i)= q'*Q(:,i);
    q = q - h(i)*Q(:,i);
  end
  h(k+1) = norm(q);
  q = q / h(k+1);
end