%  Use geometric multigrid to solve 1D (or 2D) Poisson equation


z = 0:N; z=z'/N; % Points on [0,1]


dx = z(2)-z(1);

n  = N-1;
e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n);


%%% Smoother
S = diag(A); S=(1.+0*S)./S; S=diag(S); sigma = 0.6666;

ue = rand(n,1);  % Exact solution (Endpoints assumed to be 0)
b = A*ue;        % RHS
x = 0*b;

ue = [0;ue; 0];

for k=1:200;

   r = b-A*x;
   for j=1:3

%      Diagnositics
       xe = [0; x; 0];
       ee = ue - xe;
       plot(z,xe,'bo-',z,ee,'r-',z,ue,'k--'); 
       legend('numeric','error','exact')
       %pause(0.1)

       dx = sigma*(S*r);
       x  = x + dx;
       r  = r - A*dx;
   end;

%  pause

end;

