
function A=exponential(n)
  %rand('seed',0);
  alpha=10^(-1/11);
  sv=1:n;
  for i=2:n
    sv(i)=alpha^(i-1);
  end
  A = orth(rand(n))*diag(sv)*orth(rand(n));
