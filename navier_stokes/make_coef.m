function[Ge,Je,Me,RX] = make_coef(Xe,Ye)

%
%  Make geometry specific metrics, etc.
%
N1=size(Xe,1); N=N1-1;
E =size(Xe,3);

[Ah,Bh,Ch,Dh,z,w] =  semhat(N);
Be = w*w';

RX=zeros(N1,N1,4,E); Je=zeros(N1,N1,E); Ge=zeros(N1,N1,3,E);Me=Je;
for e=1:E; 
   X=Xe(:,:,e); Y=Ye(:,:,e); 

   Xr  = Dh*X; Xs = X*Dh';
   Yr  = Dh*Y; Ys = Y*Dh';
   Jac = Xr.*Ys - Xs.*Yr;
   rx  =  Ys./Jac; ry  = -Xs./Jac;
   sx  = -Yr./Jac; sy  =  Xr./Jac;

   Grr = Jac.*Be.*(rx.*rx + ry.*ry);
   Grs = Jac.*Be.*(rx.*sx + ry.*sy);
   Gss = Jac.*Be.*(sx.*sx + sy.*sy);

   Je(:,:,e) = Jac;
   Me(:,:,e) = Jac.*Be;
   Ge(:,:,1,e) = Grr; Ge(:,:,2,e) = Grs; Ge(:,:,3,e) = Gss;
   RX(:,:,1,e) = rx;  RX(:,:,2,e) = sx;  RX(:,:,3,e) = ry; RX(:,:,4,e) = sy;

end;

