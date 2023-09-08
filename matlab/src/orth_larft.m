%
   function [ t ] = orth_larft( T, V, tau )
%
      m = size(V,1);
      j = size(V,2);
%
      t = zeros(j-1,1);
%
      t(1:j-1,1) = V(1:m,1:j-1)' * V(1:m,j);
      t(1:j-1,1) = -t(1:j-1,1) * tau;
      t(1:j-1,1) = T(1:j-1,1:j-1) * t(1:j-1,1);
%
   end

