function A=break9(n)
  e=1:n;
  e(end-9:end)=10^(-9);
  A=orth(rand(n))*diag(e)*orth(rand(n));
