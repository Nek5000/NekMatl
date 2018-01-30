function[ux,uy] = gradl(uL,Dh,RX)

E=size(uL,3); N1=size(uL,1); N=N1-1;

ux=zeros(N1,N1,E); uy=ux;
for e=1:E; 

   U=uL(:,:,e);

   Ur  = Dh*U; Us = U*Dh';

   ux(:,:,e) = RX(:,:,1,e).*Ur + RX(:,:,2,e).*Us;
   uy(:,:,e) = RX(:,:,3,e).*Ur + RX(:,:,4,e).*Us;

end;

