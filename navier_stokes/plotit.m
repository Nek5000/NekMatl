function  plotit(U,X,Y,Title)  % Contour plot

persistent Xf   % X on fine mesh
persistent Yf   % Y on fine mesh
persistent Jh   % Interpolation from coarse to fine

N1=size(U,1); E=size(U,3);

if isempty(Jh);
   N = N1-1; Nf = 3*N; [zf,wf]=zwgll(Nf); [z,w]=zwgll(N);
   Jh = interp_mat(zf,z);
   N1=Nf+1;
   Xf=zeros(N1,N1,E); Yf=zeros(N1,N1,E);
   for e=1:E;
      Xf(:,:,e) = Jh*X(:,:,e)*Jh';
      Yf(:,:,e) = Jh*Y(:,:,e)*Jh';
   end;
end;

for e=1:E;
   Uf(:,:,e) = Jh*U(:,:,e)*Jh';
end;


umn=glmin(Uf);  umx=glmax(Uf); du=(umx-umn)/35;
xmn=glmin(Xf);  xmx=glmax(Xf);
ymn=glmin(Yf);  ymx=glmax(Yf);

hold off
for e=1:E; 
%  mesh(Xf(:,:,e),Yf(:,:,e),Uf(:,:,e));
   contour(Xf(:,:,e),Yf(:,:,e),Uf(:,:,e),[umn:du:umx]); 
   axis([xmn xmx ymn ymx]); axis square; hold on
end;
%  title(Title)

