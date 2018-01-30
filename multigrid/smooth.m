

%  Illustration of Damped-Jacobi smoothing effect on error

z = 0:N; z=z'/N; h = z(2)-z(1); % Points and gridspacing on [0,1]

n  = N-1; e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n) / (h*h);
D = diag(diag(A)); sigma=0.66666;    % Smoother

ex= rand(n,1);   % Exact solution (Endpoints assumed to be 0)
f = A*ex;        % RHS = A*exact
ex= [0; ex; 0];  % Extend xxact solution to boundaries (for plotting)

u = 0*f;         % Initial guess
for j=1:200% Apply smoother 3 times

   u = u + sigma*(D\(f-A*u));  % Smoothing step

%  Diagnositics: Use known solution to show error behavior.
   ue = [0; u; 0]; ee =ex-ue; plot(z,ue,'bo-',z,ee,'r-',z,ex,'k--',z,0*ex,'k-'); 
   xlabel('x'); ylabel('Solution / Error');
   legend('numeric','error','exact'); norm(ee); pause(0.1)
   axis square; %xlabel('- x -'); %ylabel('Solution / Error');
   if j==3;   title('Error after 3 smoothings'); figure; end;
   if j==10;  title('Error after 10 smoothings'); figure; end;
   if j==200; title('Error after 200 smoothings'); figure; end;
end;


