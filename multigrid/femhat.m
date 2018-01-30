      function[Ah,Bh,Ch,Dh,z,w] =  femhat(N)
%
%     1D fem Stiffness, Mass, and Convection matrix analogies to SEMhat
%

      dz = 2/N; dzi=1/dz;
      z  = -1:dz:1;  z=z';

      n  = N+1;

      Ep = diag(ones(n-1,1),1);  Ep = sparse(Ep); 
      Em = Ep';
      D  = speye(n);

      Dh    = dzi*(Ep-D);
      Bh    = (Ep+Em)/6 + 2*D/3; Bh(1,1)=1/3; Bh(n,n)=1/3; Bh=dz*Bh;
      w     = Bh*ones(n,1);

      Ah    = dzi*(2*D-Ep-Em); Ah(1,1)=dzi; Ah(n,n)=dzi;
      Ch    = Bh*Dh;

      w     = Bh*ones(n);
