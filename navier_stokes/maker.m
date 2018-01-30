function[R] =  maker(Q,E,N)

%
%  Make the restriction matrix
%
%  Assumes periodicity in x direction
%
%  Numbers top and bottom rows last to make Dirichlet BCs easy
%
%
%        +--- numbered 2nd to last
%        |
%        v
%    +---------+---------+---------+---------+---------P
%    |         |         |         |         |         P <-- peiodic
%    |         |         |         |         |         P
%    |  e=1    |  e=2    |  e=3    |  ...    |  e=E    P
%    |         |         |         |         |         P
%    +---------+---------+---------+---------+---------P
%        ^
%        |
%        +--- numbered last
%


[nL,nb]=size(Q); n_dirichlet = 2*E*N;
R=speye(nb); R=R(1:(nb-n_dirichlet),:);

