dT = conv(T,ML,RX,Dh,u,v,Q); % dT = CT  (with mass matrix)

N1=size(RX,1); N=N1-1; E=size(RX,4);

T=Q*T; T=reshape(T,N1,N1,E);

for e=1:E







