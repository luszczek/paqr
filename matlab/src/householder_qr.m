%
    function[ V, R, T ] = householder_qr( A, orth_qr2 );
%
        [m,n] = size(A);
%
        V = zeros(m,n);
        R = zeros(n,n);
        T = zeros(n,n);
%
        for i = 1:min(m,n),
%
            %   construct V(1:m,i), R(i,i), T(i,i)
            [ V(i:m,i), R(i,i), T(i,i) ] = orth_qr2( A(i:m,i) );

            %   set upper-part of R(1:i-1,i)
            R(1:i-1,i) = A(1:i-1,i);

            %   construct T(1:i-1,i)
            [ T(1:i-1,i) ] = orth_larft( T(1:i-1,1:i-1), V(1:m,1:i), T(i,i) );

            %   update trailing submatrix A(i:m,i+1:n)
            [ A(i:m,i+1:n) ] = orth_larfb( A(i:m,i+1:n), V(i:m,i), T(i,i) );
%
        end
%
