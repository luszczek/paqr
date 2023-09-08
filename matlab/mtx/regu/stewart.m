function A=stewart(n);

  smax = 1;
  smin = 1e-3;
  nhalf = floor(n/2);
  v = zeros(n,1);
  v(1:nhalf) = fliplr(logspace( log10(smin), log10(smax), nhalf) ); 
  A = orth(rand(n)) * diag(v) * orth(rand(n)) + 0.1 * smin * rand(n);
  
