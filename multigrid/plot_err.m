function[ip] = plot_err(u,ex,z,V,L,str,v,ihold,ispec); 

ee=ex-u; sp=V'*ee;  %% Spectrum of error

norm(ee)

ee=[0; ee; 0]; % Prolongate error to boundary, for plotting.

hold off; if ihold==1; hold on; end; 

if ispec==0; 
       plot(z,ee,str,z,0*z,'k-'); 
       strg=sprintf('%s','Error after ',int2str(v),' V-cycles');
       xlabel('x'); ylabel('Solution / Error'); title(strg); axis square
else;
       plot(L,sp,str,L,0*L,'k-'); 
end;

ip=1;

if ihold==1; pause; end;

