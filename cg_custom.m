function [x,iters,e_res] = cg_custom(A,b,tol)
% CONJGRAD  Conjugate Gradient Method.
%   X = CONJGRAD(A,B) attemps to solve the system of linear equations A*X=B
%   for X. The N-by-N coefficient matrix A must be symmetric and the right
%   hand side column vector B must have length N.
%
%   X = CONJGRAD(A,B,TOL) specifies the tolerance of the method. The
%   default is 1e-10.
%
% Example (highlight lines between %{ and %}, then press F9):
%{
  n = 6000;
  m = 8000;
  A = randn(n,m);
  A = A * A';
  b = randn(n,1);
  tic, x = conjgrad(A,b); toc
  norm(A*x-b)
%}
% By Yi Cao at Cranfield University, 18 December 2008
% Updated on 6 Feb 2014.
%   
    e_res = [];
    x_ans = A\b;
    norm_x = (x_ans'*A*x_ans)^0.5;
    if nargin<3
        tol=1e-10;
    end
    x = zeros(size(b));
    e_res = [e_res ((x_ans-x)'*A*(x_ans-x))^0.5/norm_x ];
    r = b - A*x;

    y = -r;
    z = A*y;
    s = y'*z;
    t = (r'*y)/s;
    x = x + t*y;
    iters = 0;
    for k = 1:numel(b)
       r = r - t*z;
       if( norm(r) < tol )
            break;
       end
       B = (r'*z)/s;
       y = -r + B*y;
       z = A*y;
       s = y'*z;
       t = (r'*y)/s;
       x = x + t*y;
       e_res = [e_res ((x_ans-x)'*A*(x_ans-x))^0.5/norm_x ];
       iters=iters+1;
    end
    
 end