      function[l] =  legendre(x,N) 

%     Compute Legendre polynomials up to deg. N at points x

      m=length(x);
      l=ones(m,N+1);
      l(:,2)=x;

      for k=2:N; i=k+1;
         l(:,i) = ( (2*k-1)*x.*l(:,i-1)-(k-1)*l(:,i-2) )/k;
      end;
