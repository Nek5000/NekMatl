
N=1;E=100;

Q=makeq(E,N);
R=maker(Q,E,N);
[X,Y]=make_geom(E,N);
[G,J,RX]=make_coef(X,Y);
[Ah,Bh,Ch,Dh,z,w]=semhat(N);

[n,nb]=size(R); n=nb;

A=zeros(n,n);
for j=1:n;
    x=zeros(n,1); x(j)=1; A(:,j)=abar(x,Dh,G,Q); 
end;

[V,D]=eig(A); d=sort(diag(D)); s=sqrt(d);
