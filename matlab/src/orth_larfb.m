%
   function [ A ] = orth_larfb( A, v, tau )
%
      [m,n] = size(A);
%
      alpha = v(1:m,1)' * A(1:m,1:n);
      alpha = tau * alpha;
      A(1:m,1:n) = A(1:m,1:n) - v(1:m,1) * alpha;
%
   end

