%
    function[ V, R, T, dead_cols] = householder_poqr( A, orth_qr2, is_deficient, nrm, nrmA );
%    function[ V, R, T, dead_cols , cri, dd] = householder_poqr( A, orth_qr2, is_deficient, nrm, nrmA );
%
        [m,n] = size(A);
        Aorig = A;
%
        dead_cols = zeros(n,1);
        V = zeros(m,n);
        R = zeros(n,n);
        T = zeros(n,n);
%cri = zeros(n,1);
%dd = zeros(n,1);
%
        k = 1;
        for i = 1:min(m,n),
%
            %   construct V(1:m,k), R(k,k), T(k,k) based on A(k:m,i)
            [ V(k:m,k), R(k,k), T(k,k) ] = orth_qr2( A(k:m,i) );
%
%cri(i) = 0.1 * eps * max(vecnorm(Aorig(1:end,1:k)));
%dd(i) = abs( R(k,k) );
           if ( is_deficient(R, k, nrm, Aorig, nrmA) )
           %if (dd(i) < cri(i))
           %if ( abs( R(k,k) ) < tol ),
                dead_cols(i) = 1;
            else
                %   set upper-part of R(1:k-1,k)
                R(1:k-1,k) = A(1:k-1,i);
%
                %   construct T(1:k-1,k)
                [ T(1:k-1,k) ] = orth_larft( T(1:k-1,1:k-1), V(1:m,1:k), T(k,k) );
%
                %   update trailing submatrix A(k:m,k+1:n)
                [ A(k:m,i+1:n) ] = orth_larfb( A(k:m,i+1:n), V(k:m,k), T(k,k) );
%
                k = k + 1;
            end
%
        end
%
        V = V(1:m,1:k-1);
%        R = R(1:k-1,1:k-1);
        R = R(1:end,1:k-1);
        T = T(1:k-1,1:k-1);
