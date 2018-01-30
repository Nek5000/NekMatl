%
%   This tests the solution to Laplaces equation for an E x 1 array
%   of spectral elements of order N.
%
%   Geometry is determined by make_geom.
%
%   Periodicity in the x-direction established by Q=makeq()
%
%   Dirichlet conditions at top and bottom established by R=maker()
%

N=10;E=9; N1=N+1; nL=N1*N1*E;

Q=makeq(E,N);
R=maker(Q,E,N);
[X,Y]=make_geom(E,N);        % Geometry in local form
xM=glmax(X); xm=glmin(X); s=2*pi/(xM-xm); X=s*X;
[G,J,RX]=make_coef(X,Y);
[Ah,Bh,Ch,Dh,z,w]=semhat(N);

[n,nb]=size(R);

UL=sin(X).*exp(Y);  


%   maskL := Convenient local (diagonal) mask to set BCs
%  amaskL := "anti-mask" -- zeros out interior values

mask=diag(R'*R); maskL=reshape(full(Q*mask),N1,N1,E); amaskL=1-maskL;
Ub=amaskL.*UL;
plotit(Ub,X,Y,'Ubar');pause(1)

FL=-aloc(Ub,Dh,G);      
f=R*Q'*reshape(FL,nL,1);  % RHS is:   f = - R*Q'*A_L*Ub_L

%  Generate Abar
Ab=zeros(nb,nb);
for j=1:nb;
    x=zeros(nb,1); x(j)=1; Ab(:,j)=abar(x,Dh,G,Q); 
end;

A=R*Ab*R';
u0 = A\f;
U0 = Q*(R'*u0); U0=reshape(U0,N1,N1,E);
U=U0+Ub;

scale=glmax(abs(U-UL))
plotit((U-UL)/scale,X,Y)

