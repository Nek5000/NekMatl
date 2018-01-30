function [uL,vL,pr]=pressure_project(ut,vt,Ab,Q,ML,RX,Dh); % Div-free velocity

persistent Ai;
if isempty(Ai); Ai=pinv(full(Ab)); end;

N1=size(ut,1); E=size(ut,3); nL=N1*N1*E;
re=ut;
for e=1:E;
  rr=ML(:,:,e).*(RX(:,:,1,e).*ut(:,:,e)+RX(:,:,3,e).*vt(:,:,e));
  rs=ML(:,:,e).*(RX(:,:,2,e).*ut(:,:,e)+RX(:,:,4,e).*vt(:,:,e));
  re(:,:,e)=Dh'*rr+rs*Dh;
end;

re=ut;
for e=1:E;
  rr=ML(:,:,e).*(RX(:,:,1,e).*ut(:,:,e)+RX(:,:,3,e).*vt(:,:,e));
  rs=ML(:,:,e).*(RX(:,:,2,e).*ut(:,:,e)+RX(:,:,4,e).*vt(:,:,e));
  re(:,:,e)=Dh'*rr+rs*Dh;
end;


r=Q'*reshape(re,nL,1); n=size(r,1); r=r-sum(r)/n;
pr=Ai*r; pr=Q*pr; pr=reshape(pr,N1,N1,E);
[prx,pry]=gradl(pr,Dh,RX);
uL = ut - prx; vL = vt - pry;

