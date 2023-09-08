% LAPACK?

function [ v, r, t ] = orth_geqr2( x )

    m = size(x,1);
    normx = norm( x(2:m,1), 2);
    norma = sqrt( normx*normx + x(1,1)*x(1,1) );
    if ( x(1,1) > 0.0e+00 )
        r = - norma; else r = norma;
    end
    x(1,1) = x(1,1) - r;
    v(1,1) = 1.0e+00;
    v(2:m,1) = x(2:m,1) / x(1,1);
    t = 2.0e+00 / ( 1.0e+00 + ( normx / x(1,1) )^2 );

    end

