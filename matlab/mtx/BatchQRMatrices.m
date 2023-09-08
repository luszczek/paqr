

% Create a matrix for Pth-order interpolation in 3D
P = 6; % 2-6 are good values to vary the size from 35x10 to 165x84
N = P+3; % extra points for LS interpolation

% Create the matrix of 3D interpolation powers
col = 1;
clear Powers;
for i=0:P
for j=0:P
for k=0:P
    if (i+j+k) <= P
        Powers(1:3,col) = [i;j;k];
        col = col+1;
    end
end
end
end
Ncol = size(Powers,2);
assert(Ncol == (P+1)*(P+2)*(P+3)/6);

% Use these powers to create the interpolation matrix
Nmat = 100;
batchA = cell(Nmat,1);
for mat=1:Nmat

row=1;
clear A;
    
% Random centering for the interpolation for each matrix
x0 = rand();
y0 = rand();
z0 = rand();
for x=1:N
    xpow = (x - x0).^Powers(1,:);
for y=1:N
    ypow = (y - y0).^Powers(2,:);
for z=1:N
    if (x+y+z) <= N+2
        zpow = (z - z0).^Powers(3,:);        
        A(row,1:Ncol) = xpow .* ypow .* zpow;
        row = row+1;
    end
end
end
end

batchA{mat,1} = A;
end

%% Nrow = size(A,1);
%% assert(Nrow == N*(N+1)*(N+2)/6);
%
%% Check if the iterative algorithm works for each one
%for mat=1:Nmat
%    mat
%    A = batchA{mat,1};
%
%    % initial guess
%    scale = norm(A,inf)*norm(A,1);
%    V = pinv(A) + A' / scale;
%    maxiter = 100;
%    m = size(A,1);
%    n = size(A,2);
%    tol = 1e-11*sqrt(scale);
%    clear err;
%    for iter = 1:maxiter
%        X = A*V; 
%        V = V*(2*eye(m,m) - X); % 2nd-order method
%%         V = V*(3*eye(m,m) - X*(3*eye(m,m) - X)); % 3rd-order method
%        % 7th-order method
%%         V = V*(7*eye(m,m) + X*(-21*eye(m,m) + X*(35*eye(m,m) + ...
%%             X*(-35*eye(m,m) + X*(21*eye(m,m) + X*(-7*eye(m,m) + X))))));
%        err(iter) = norm(A*V*A - A,'fro');
%        if (err(iter) < tol)
%            break;
%        end
%    end
%    iter
%    err
%end

