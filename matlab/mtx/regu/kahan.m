function A=kahan(n)

  c=0.6;
  s=sqrt(1-0.36);
  A=ones(n);
  B=ones(n);
  for i=2:n
    A(i,i)=s^(i-1);
    B(i,i+1:end)=-c;
  end
  A=A*B;

