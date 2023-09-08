% like DLARFG in LAPACK

function [ v, r, t ] = house_lapack (x)
    v = x;
    n = length (x) ;
    r = norm (x) ;
    if (x (1) > 0)
        r = -r ;
    end
    t = (r - x (1)) / r ;
    v (2:n) = x (2:n) / (v (1) - r) ;
    v (1) = 1 ;

