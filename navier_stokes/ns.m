
format compact; format longe; %figure


[X,Y,G,RX,Jac,Area,Nhat]=mygeom(E,N);
[Ah,Bh,Ch,Dh,zpts,wt]=semhat(N);

dxmin=mindx(X,Y);   % Get min dx for CFL constraint

Tfinal=1; dt=CFL*dxmin; nstep=ceil(Tfinal/dt); dt=Tfinal/nstep;

for step=1:nstep; time=step*dt;

    if step==1; b0=1.0;    b= [ -1 0 0 ]';       a=[ 1  0 0 ]'; end;
    if step==2; b0=1.5;    b=([ -4 1 0 ]')./2;   a=[ 2 -1 0 ]'; end;
    if step==3; b0=11./6.; b=([ -18 9 -2 ]')./6; a=[ 3 -3 1 ]'; end;
    b0i = 1./b0;

%   Nonlinear step

    fx1 = -conv(u,J,bo,dh,uo,uo,Q,E); % du = Cu  (with mass matrix)
    fy1 = -conv(v,J,bo,dh,uo,uo,Q,E); % dv = Cv
    fy1 = fy1 + BL(T,J);              % Thermal forcing (Boussinesq)

    rx  = a(1)*fx1+a(2)*fx2+a(3)*fx3; % kth-order extrapolation
    ry  = a(1)*fy1+a(2)*fy2+a(3)*fy3; % kth-order extrapolation
    fx3=fx2; fx2=fx1; fy3=fy2; fy2=fy1;

    rx  = dt*rx - (b(1)*u+b(2)*u2+b(3)*u3); u3=u2; u2=u; % Add BDF terms
    ry  = dt*ry - (b(1)*v+b(2)*v2+b(3)*v3); v3=v2; v2=u;

    ut  = b0i*(Bi.*rx); vt  = b0i*(Bi.*ry);  % Velocity after nonlinear step (no bc)

    [u,v,pr] = pressure_project(ut,vt,G,Jac,...,pmask); % Div-free velocity

    

format compact; format longe; %figure

No = N;
No = N+1; 

vtype = 'c ';  % Plane-rotation velocity field
vtype = '+ ';  % Straining velocity field

[ah,bh,ch,dh,zpts,wt]=semhat(N);
[ao,bo,co,do,zpto,wo]=semhat(No);
J = interp_mat(zpto,zpts);

E = 3; %  E = input ('Input E: ');

Len = 1.0;
Le  = Len/E;

nn = E*N+1;
xn = zeros(nn,1); x0 = 0; j1 = 1;
for e=1:E; j0 = j1; j1 = j0 + N;
    xn(j0:j1) = x0 + (Le/2)*(zpts+1); x0 = xn(j1);
end;
x = xn(2:j1-1); n = nn-2;

xo = zeros(E*No+1,1); x0 = 0; j1 = 1;
for e=1:E; j0 = j1; j1 = j0 + No;
    xo(j0:j1) = x0 + (Le/2)*(zpto+1); x0 = xo(j1);
end;

hold off; [vx,vy,vxq,vyq,bx,sx,by,sy] = velocity(xo,xo,vtype,'T');

Qt = buildQ(E,No); 
vxo = glob2loc(vxq,No,E,E,Qt,Qt);
vyo = glob2loc(vyq,No,E,E,Qt,Qt);
Q = buildQ(E,N); n = size(Q,2); Q = Q(:,2:n-1);  %  Assign Dirichlet BCs

n = nn-2;
Cf = sparse(zeros(n*n,n*n));
u  = zeros(n*n,1);

for j=1:n*n;
    u(j) = 1;
    Cf(:,j) = Conv(u,J,bo,dh,vxo,vyo,Q,E,1.);
    u(j) = 0;
end;

Cf=full(Cf);
[vc,lc]=eig(Cf); lc=diag(lc);

pause
plot(lc,'ro'); 
pause
axis equal
pause



nn = E*N+1;
an = zeros(nn,nn); bn = zeros(nn,nn); x0 = 0; j1 = 1;
for e=1:E; j0 = j1; j1 = j0 + N;
    an(j0:j1,j0:j1) = an(j0:j1,j0:j1) + (2/Le)*ah;
    bn(j0:j1,j0:j1) = bn(j0:j1,j0:j1) + (Le/2)*bh;
end;
a = an(2:j1-1,2:j1-1); b = bn(2:j1-1,2:j1-1); n = nn-2;

b = diag(b);
bi= 1./b;
Bi= bi*bi';
Bi= reshape(Bi,n*n,1);

eh=ones(n,1);
xq=x*eh'; yq=eh*x';
[u,uq] = uinit(x,x);

f=zeros(n*n,3);
dx = diff(x); dxmin = min(dx);
dt = CFL*dxmin; Tfinal = 1; nstep = 5*Tfinal/dt; %nstep=500;
[vx,vy,vxq,vyq,bx,sx,by,sy] = velocity(xo,xo,vtype,'T');
hold on;
uuuu
%
%   Single domain, explicit convection, no dealiasing, no filter
%
%   NOTE: w/o dealiasing, this case is restricted to
%         constant velocities or to rotational flows.
%
%   The purpose here is to demonstrate a minimalist SEM
%   implementation of convection that is as fast as possible.
%
%   For this particular case (single element), no mass matrix
%   is required.

[Ah,Bh,Ch,Dh,z,w] = SEMhat(N); 


nx=size(Ah,1); x=z; 
Ax=eye(nx); Ax(2:nx-1,2:nx-1)=Ah(2:nx-1,2:nx-1);  % Interior only, for
Bx=eye(nx); Bx(2:nx-1,2:nx-1)=Bh(2:nx-1,2:nx-1);  % Dirichlet conditions
bx=diag(Bx); Bx=sparse(Bx);
Dx=Dh;

y=x; ny=nx; Ay=Ax; By=Bx; Dy=Dx; by=bx;

% Diagonal mask (1 in interior, 0 on boundary)
mx=ones(nx,1); mx(1)=0; mx(nx)=0;
my=ones(ny,1); my(1)=0; my(ny)=0;
mask = mx*my';


%  Define Mesh
[xx,yy] = meshgrid(x,y); xx=xx';yy=yy'; % Lexicographical ordering

pi2=pi/2;    %  Define Velocity
%cx = -cos(pi2*xx).*sin(pi2*yy); cy =  sin(pi2*xx).*cos(pi2*yy); 
%cx = .5+0*yy; cy=1-0*xx; 
 cx = -yy; cy=xx; 


% Premultiply velocity by mask to enforce Dirichlet condition
 cxo = mask.*cx; cyo = mask.*cy;


% Set initial condition, zero on boundary
x0 =-0.0; y0=-0.5; delta = 0.10; rr=(xx-x0).^2+(yy-y0).^2;
u0 = exp(-((rr./(delta^2)).^1)).*(1-xx.*xx).*(1-yy.*yy);
mesh(xx,yy,u0); pause(.1);

Tfinal = 2*pi;

dxmin = min(abs(diff(x))); Dt = 0.0625*dxmin;    % CFL constraint
Nsteps = ceil(Tfinal/Dt); Dt = Tfinal/Nsteps; % Integer # steps/period


u=u0; u1=0*u; u2=u1; e1=u1; e2=u1; 


Dt
Nsteps

time=0;  Nstep=20;
for step=1:Nstep; time=step*dt;

  if step==0; 
uuuu
uuuu

    %Save O(dt) soln at dt and 2*dt:
    u11=u1; u01=u; e21=e2; e11=e1; 

    istart=1; istop=4; time = 0; dt=Dt/2.;
    u=u0; u1=0*u; u2=u1; e1=u1; e2=u1; 
    b0=1; b=[ -1 0 0 ]'; a=[ 1 0 0 ]';
  elseif step==1; 

    %Save O(dt/2) soln at dt and 2*dt:
    u12=u1; u02=u; e22=e2; e12=e1; 

    % Richardson Ext.
    u2=u0; u1=(4*u12-u11)/3; u=(4*u02-u01)/3;
    e2=(4*e22-e21)/3; e1=(4*e12-e11)/3;

    b0=11./6.; b=[ -18 9 -2 ]; b=b'./6; a=[ 3 -3 1 ]';

  end;

  b0i = 1./b0;

  for step=istart:istop; time=time+dt;

    ux = cxo.*(Dx*u); uy = cyo.*(u*Dy');
    e3=e2; e2=e1; e1=ux+uy;

    u3=u2; u2=u1; u1=u;

    u = -b0i*( u1*b(1)+u2*b(2)+u3*b(3) + dt*(e3*a(3)+e2*a(2)+e1*a(1)));

%   u = F*u*F';                       %  FILTER high wavenumbers
        
%   if mod(step,100)==0; mesh(xx,yy,u); pause(.1); end;
    if mod(step,1000)==0; mesh(xx,yy,u); axis([-1 1 -1 1 -1 1]); pause(.1); end;
    if mod(step,Nsteps)==0; err = u-u0; max(max(abs(err))), end;

  end;
  mesh(xx,yy,u); axis([-1 1 -1 1 -1 1]); pause(.1);
end;
err = u-u0; max(max(abs(err)))
%mesh(xx,yy,err); pause(.1);
uuuu
%
%   Convect a Gaussian pulse using plane rotation on [-1,1]^2.
%
%   Usage:   
%
%   N=33; cfast
%
%   Single domain, explicit convection, no dealiasing, no filter
%
%   NOTE: w/o dealiasing, this case is restricted to
%         constant velocities or to rotational flows.
%
%   The purpose here is to demonstrate a minimalist SEM
%   implementation of convection that is as fast as possible.
%
%   For this particular case (single element), no mass matrix is required.
%
%   This version has only an O(dt^2) accurate initial step.
%   See c_acc.m for a more accurate startup phase.
%
%

[Ah,Bh,Ch,Dh,z,w] = SEMhat(N); 

x=z; Dx=Dh; nx=size(Ah,1); 
y=z; Dy=Dh; ny=size(Ah,1); 


% Diagonal mask (1 in interior, 0 on boundary)
mx=ones(nx,1); mx(1)=0; mx(nx)=0;
my=ones(ny,1); my(1)=0; my(ny)=0;
mask = mx*my';


%  Define Mesh
[xx,yy] = meshgrid(x,y); xx=xx';yy=yy'; % Lexicographical ordering

pi2=pi/2;    %  Define Velocity  (comment/uncomment to test other velocities)
%cx = -cos(pi2*xx).*sin(pi2*yy); cy =  sin(pi2*xx).*cos(pi2*yy); 
%cx = .5+0*yy;                   cy=1-0*xx; 
 cx = -yy;                       cy=xx; 


% Premultiply velocity by mask to enforce Dirichlet condition
 cxo = mask.*cx; cyo = mask.*cy;


% Set initial condition; enforce zero on boundary w/ quadratic bubble
x0 =-0.0; y0=-0.5; delta = 0.10; rr=(xx-x0).^2+(yy-y0).^2;
u0 = exp(-((rr./(delta^2)).^1)).*(1-xx.*xx).*(1-yy.*yy);
mesh(xx,yy,u0); axis([-1 1 -1 1 -.1 1]); pause(.1);

Tfinal = 2*pi;                                   % Time for single revolution

CFL = .5;                                        % CFL constraint ( < ~.5 )
dxmin = min(abs(diff(x))); Dt = CFL*dxmin;    
Nsteps = ceil(Tfinal/Dt);  Dt = Tfinal/Nsteps;   % Integer # steps/period
format long; [Dt Nsteps]


Npass=20; time=0; dt=Dt;

u=u0; u1=0*u; u2=u1; e1=u1; e2=u1;   % Initialize lagged fields

for pass=1:Npass;                    % Simple scheme, only
  for step=1:Nsteps;                 %    O(dt^2) accuracy at start



    ux = cxo.*(Dx*u); uy = cyo.*(u*Dy');
    e3=e2; e2=e1; e1=ux+uy;

    u3=u2; u2=u1; u1=u;

    u = -b0i*( u1*b(1)+u2*b(2)+u3*b(3) + dt*(e3*a(3)+e2*a(2)+e1*a(1)));

  end;

  mesh(xx,yy,u); axis([-1 1 -1 1 -.1 1]); pause(.1);
  err(pass) = max(max(abs(u-u0))); t(pass) = time; 
  [pass time err(pass)]

end;

figure
semilogy(t,err,'ro',t,err,'r-'); axis square;
xlabel('time'); ylabel('maximum pointwise error');
title('Error for Convection, Std. Start')

