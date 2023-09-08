% Higham

function [ v, r, t ] = house_higham (x)

    n = length (x) ;
    if (n == 1)
        sigma = 0 ;
    else
        sigma = x (2:n)' * x (2:n) ;
    end
    v = x ;
    if (sigma == 0)
        r = x (1) ;
        v (1) = 0 ;
        t = 0 ;
    else
        r = sqrt (x(1)^2 + sigma) ;
        if (x (1) <= 0)
            v (1) = x (1) - r ;
        else
            v (1) = -sigma / (x (1) + r) ;
        end
        t = -1 / (r * v(1)) ;
    end

