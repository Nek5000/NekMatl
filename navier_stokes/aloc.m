function[wL] = abar(uL,Dh,G)

E=size(G,4); N1=size(G,1); N=N1-1;

wL=uL;

for e=1:E; 

   ur = Dh*uL(:,:,e); us = uL(:,:,e)*Dh';
   wr = G(:,:,1,e).*ur + G(:,:,2,e).*us;
   ws = G(:,:,2,e).*ur + G(:,:,3,e).*us;
   wL(:,:,e) = Dh'*wr + ws*Dh;

end;

