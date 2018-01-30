function[X,Y] =  make_geom(E,N)

[z,w]=zwgll(N);

N1=N+1;

[X0,Y0]=ndgrid(z,z); X=zeros(N1,N1,E); Y=X;

x0=0; x1=x0+1; y0=0; y1=y0+1;
for e=1:E;  % Domain of unit [0,1] squares
    X(:,:,e) = x0+.5*(x1-x0)*X0; x0=x0+1; x1=x0+1;
    Y(:,:,e) = y0+.5*(y1-y0)*Y0; 
end;


xM=glmax(X); xm=glmin(X); s=2*pi/(xM-xm); X=s*X; %  X on [0,2pi]


%  Map to circle, starting at theta=pi/2
xmax=glmax(X); xmin=glmin(X);
ymax=glmax(Y); ymin=glmin(Y);
theta = (pi/2)-(2*pi)*(X-xmin)/(xmax-xmin);
rad   = 0.5 + (Y-ymin)/(ymax-ymin);
X=rad.*cos(theta); Y=rad.*sin(theta);

%% 
%% Below just for diagnostics
%% 
%% F=X.*X-Y.*Y; fmx=max(max(max(F))); fmn=min(min(min(F)));
%% 
%% hold off
%% for e=1:E; 
%%    mesh(X(:,:,e),Y(:,:,e),F(:,:,e)); hold on
%%    contour(X(:,:,e),Y(:,:,e),F(:,:,e),[fmn:.1:fmx]); hold on
%% end;

