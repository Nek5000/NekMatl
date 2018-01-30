
N=10;E=5; N1=N+1; nL=N1*N1*E;

Q=makeq(E,N);
R=maker(Q,E,N);
[X,Y]=make_geom(E,N);        % Geometry in local form
xM=glmax(X); xm=glmin(X); s=2*pi/(xM-xm); X=s*X;
[G,J,ML,RX]=make_coef(X,Y);
[Ah,Bh,Ch,Dh,z,w]=semhat(N);

[n,nb]=size(R);


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

uuuu

dxmin=mindx(X,Y);   % Get min dx for CFL constraint
Tfinal=1; dt=CFL*dxmin; nstep=ceil(Tfinal/dt); dt=Tfinal/nstep;

for step=1:nstep; time=step*dt;

    if step==1; b0=1.0;    b= [ -1 0 0 ]';       a=[ 1  0 0 ]'; end;
    if step==2; b0=1.5;    b=([ -4 1 0 ]')./2;   a=[ 2 -1 0 ]'; end;
    if step==3; b0=11./6.; b=([ -18 9 -2 ]')./6; a=[ 3 -3 1 ]'; end;
    b0i = 1./b0;

%   Nonlinear step - unassembled, not multiplied by mass matrix

    fx1 = -conv(u,RX,Dh,u,v); % du = Cu  (without mass matrix)
    fy1 = -conv(v,RX,Dh,u,v); % dv = Cv  (without mass matrix)
    fy1 = fy1 + T;            % Thermal forcing (Boussinesq)
    ft1 = -conv(T,RX,Dh,u,v); % dT = CT  (without mass matrix)

    rx  = a(1)*fx1+a(2)*fx2+a(3)*fx3; % kth-order extrapolation
    ry  = a(1)*fy1+a(2)*fy2+a(3)*fy3;
    rt  = a(1)*ft1+a(2)*ft2+a(3)*ft3;
    fx3=fx2; fx2=fx1; fy3=fy2; fy2=fy1; ft3=ft2; ft2=ft1;

    rx  = dt*rx - (b(1)*u+b(2)*u2+b(3)*u3); u3=u2; u2=u; % Add BDF terms
    ry  = dt*ry - (b(1)*v+b(2)*v2+b(3)*v3); v3=v2; v2=v;
    rt  = dt*rt - (b(1)*T+b(2)*T2+b(3)*T3); T3=T2; T2=T;

    ut  = b0i*rx; vt = b0i*ry; T = b0i*rt; % Tentative fields

    [uL,vL,pr]=pressure_project(ut,vt,G,Jac,...,pmask); % Div-free velocity

    u = R*(Q'*(ML*uL)); u = UH\(LH\(R*u); u = Q*(R'*u); u=reshape(u,N1,N1,E);
    v = R*(Q'*(ML*vL)); v = UH\(LH\(R*v); v = Q*(R'*v); v=reshape(v,N1,N1,E);
    T = R*(Q'*(ML*T));  T = UH\(LH\(R*T); T = Q*(R'*T); T=reshape(T,N1,N1,E);

end;
