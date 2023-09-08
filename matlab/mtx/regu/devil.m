function A=devil(n)
  l=20;
  s=zeros(n,1);
  nst=floor(n/l);
  for i=1:nst
    s(1+l*(i-1):l*i)=-0.6*(i-1);
  end
  s(l*nst:end)=-0.6*(nst-1);
  s=10.^s;
  A=orth(rand(n))*diag(s)*orth(randn(n));
