%%laplace.m
%% 2-d laplace-poisson equation, set up as in Roache, p. 113
%% for a rectangular domain -M grid points
%% in x direction, N in y-direction
%% Deltax=Deltay. Setting Deltax determines the physical dimensions
%%C Wunsch, 1989. Modified 1999

clear
N = 5
M = 5
Deltax = 1
Deltay = 1

A=zeros(N,M);

%set up matrix A for interior values first.
%%numbering goes along rows starting at top edge:
    for j9=2:N-1
   for i9=2:M-1
k9=(j9-1)*M+i9; %this is the grid point number
%%find neighbors:
nbhs=[k9-1,k9+1,k9-M,k9+M]; %left, right, above, below
A(k9,nbhs)=ones(1,4);A(k9,k9)=-4;
   end; %end i9 loop
   end; %end j9 loop

%%now set equations of boundary conditions
%top boundary:
A(1:M,1:M)=eye(M);
%bottom boundary:
A(N*M-M+1:N*M,N*M-M+1:N*M)=eye(M);
%left boundary:
A(M+1:M:(N*M-2*M+1),M+1:M:(N*M-2*M+1))=eye(N-2);
%right boundary:
A(2*M:M:(N*M-M),2*M:M:(N*M-M))=eye(N-2);

%%now set up right hand-side: bc's
% are set and interior source values are set:
b=sparse(zeros(M*N,1));
b(1:M,1)=zeros(M,1);%top bcs; change to suit problem
b(N*M-M+1:N*M,1)=zeros(M,1); %bottom bc
b(M+1:M:(N*M-2*M+1),1)=zeros(N-2,1); %left boundary conds.
b(2*M:M:(N*M-M),1)=zeros(N-2,1); %right bcs:
%%interior sources: change to suit problem
for j9=2:N-1
  for i9=2:M-1
  k9=(j9-1)*M+i9; %interior grid no. again
  b(k9)=0*(Deltax)^2; %%note appearance of Deltax. All 0 here.
   end, %end i9 loop
   end, %end j9 loop

b(13)=1*Deltax^2; %special case

