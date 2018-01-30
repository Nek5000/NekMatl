function[ax] = abar(u,Dh,G,Q)

E=size(G,4); N1=size(G,1); N=N1-1;

ax=u; uL=reshape(Q*u,N1,N1,E); wL=uL;

for e=1:E; 

   ur = Dh*uL(:,:,e); us = uL(:,:,e)*Dh';
   wr = G(:,:,1,e).*ur + G(:,:,2,e).*us;
   ws = G(:,:,2,e).*ur + G(:,:,3,e).*us;
   wL(:,:,e) = Dh'*wr + ws*Dh;

end;
wL=reshape(wL,N1*N1*E,1);
ax=Q'*wL;

