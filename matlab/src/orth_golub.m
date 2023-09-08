% Golub

function [ v, r, t ] = house_golub (x)

    m = size(x);
    if (m == 1)
        sigma = 0 ;
    else
        sigma = x (2:m)' * x (2:m) ;
    end
    v = [1; x(2:m)];
    if (sigma == 0)
        t = 0;
    else
        s = sqrt(x(1)*x(1) + sigma);
        if (x(1) <= 0)
            v(1) = x(1) - s;
        else
            v(1) = -sigma/(x(1) + s);
        end
        t = 2*v(1)*v(1)/(sigma + v(1)*v(1));
        v = v / v(1);
    end
    r = x(1) - t*v(1)*v'*x;

