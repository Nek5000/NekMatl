
N=7;E=5;

Q=makeq(E,N);
R=maker(Q,E,N);
[X,Y]=make_geom(E,N);        % Geometry in local form
[G,J,RX]=make_coef(X,Y);
[Ah,Bh,Ch,Dh,z,w]=semhat(N);

[n,nb]=size(R); n=nb;

Ab=zeros(nb,nb);
for j=1:n;
    x=zeros(n,1); x(j)=1; Ab(:,j)=abar(x,Dh,G,Q); 
end;

[V,D]=eig(Ab); d=sort(diag(D)); s=sqrt(d);

norm(Ab-Ab')


UL=exp(X).*sin(Y); 


RL=FL;
for e=1:E;
    RL(:,:,e)=J(:,:,e).*(Bh*FL(:,:,e)*Bh');
end;
r=




