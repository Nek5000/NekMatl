
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Geometric multigrid to solve 1D Poisson equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Usage:   N=40; m=1; mgrid
%
%  N = number of points-1, including domain endpoints.
%      Number of degrees of freedom is N-1.
%
%  k = number of fine-grid smoothings per V-cycle
%

 hold off; format compact;  close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Illustration of Damped-Jacobi smoothing effect on error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 z = 0:N; z=z'/N; h = z(2)-z(1); % Points and gridspacing on [0,1]
 n = N-1; e = ones(n,1);
 A = spdiags([-e 2*e -e], -1:1, n, n) / (h*h);

% E=1; Lx=1; bc=0; [A,Bsem,Csem,z] = sem1c(E,N,Lx,bc); z=[0; z ; 1];

 D = diag(diag(A)); sigma=0.66666;    % Smoother
 sigma=0.33333;    % Smoother


 Nc = N/2; zc = 0:Nc; zc=zc'/Nc;     %% Build Level 2 Operators:
 P = 0*A; i=0;                       %% 
 for j=1:(n/2); i=i+1; P(i,j)=.5; i=i+1; P(i,j)=1; P(i+1,j)=.5; end;
               i=i+1; P(i,j)=.5;

 P=P(1:i,1:j);  %% PROLONGATION OPERATOR
 R=P';          %% RESTRICTION OPERATOR
 Ac = R*A*R';   %% COARSE-GRID SYSTEM


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [V,L]=eig(full(D\A)); 
 L=eig(L);  %% Note: eigenfunctions indeterminate up to a sign
 S=V'*V; S=diag(S); S=1./sqrt(S); S=diag(S); V=V*S; 

 for k=1:n; if V(1,k) < 0; V(:,k)=-V(:,k); end; end;

 xhat=ones(n,1);
%ex = V*xhat;
 ex= rand(n,1); 


 f = A*ex;                      % f = A*exact
 u = 0*f; r = f-A*u;            % Initial guess and residual

 v  = 0; ispec=0; % ispec=1 --> plot spectra
 ip = plot_err(u,ex,z,V,L,'g-',v,0,ispec); %pause

 m=2;
 for v=1:10;                    % 10 V-cycles
   for j=1:m                    % Apply smoother m times
     u=u+sigma*(D\r); r=f-A*u;
   end;
   ip = plot_err(u,ex,z,V,L,'b-',v,1,ispec); 

   ef = R'*(Ac \ (R*r) );       % Coarse-grid correction
   u  = u + ef;
   ip = plot_err(u,ex,z,V,L,'r-',v,1,ispec);  %pause

 end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
