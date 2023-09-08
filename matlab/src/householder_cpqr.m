%
    function[ V, R, T, P ] = householder_cpqr( A, orth_qr2 );
%
        [m,n] = size(A);
%
        nrms = zeros(n,1);
        P_index = [1:n];
%
        V = zeros(m,n);
        R = zeros(n,n);
        T = zeros(n,n);
        P = eye(n,n);
%
        for i = 1:min(m,n),
%
            nrms(i:n) = vecnorm(A(i:m,i:n),2);
            [~,idx] = max( nrms(i:n) );
            idx = idx+i-1;
%
            if i < idx
%
	       tmp = P_index(i);
               P_index(i) = P_index(idx);
               P_index(idx) = tmp;
%
               tmp = A(1:m,i);
               A(1:m,i) = A(1:m,idx);
               A(1:m,idx) = tmp;
%
            end
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
        P = P( :, P_index );
%
