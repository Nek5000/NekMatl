function [dT] = convl(T,RX,Dh,u,v); % dT = CT  (with mass matrix)

N1=size(T,1); N=N1-1; E=size(T,3);

dT=T;

for e=1:E
    Tr = Dh*T(:,:,e); Ts = T(:,:,e)*Dh';
    Tx = Tr.*RX(:,:,1,e) + Ts.*RX(:,:,2,e);
    Ty = Tr.*RX(:,:,3,e) + Ts.*RX(:,:,4,e);
    dT(:,:,e) = (u(:,:,e).*Tx + v(:,:,e).*Ty);
end;
