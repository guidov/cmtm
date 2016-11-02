
%%Neumann.m
%% 2-d laplace-poisson equation, set up as in Roache, p. 113
%%with Neumann boundary conditions. Cf Laplace.m
%% for a rectangular domain -M in x direction, N in y-direction
%% M is number of grid points in EACH direction. (1,1) is northwest
%%corner, N,1 is sw. corner, etc. Deltax=Deltay
%%C Wunsch 1989. modified 1999
A=zeros(N*M);
%set up matrix A for interior values first.
%%numbering goes along rows starting at top edge:
count=0;intind=zeros(1,(N-2)*(M-2));
    for j9=2:N-1
   for i9=2:M-1
k9=(j9-1)*M+i9; %this is the grid point number
count=count+1;
intind(1,count)=k9;  %interior index
%%find neighbors:
nbhs=[k9-1,k9+1,k9-M,k9+M]; %left, right, above, below
A(k9,nbhs)=ones(1,4);A(k9,k9)=-4;
   end; %end i9 loop
   end; %end j9 loop
%%now set equations of boundary conditions; + OUTWARD
%top boundary:
A(2:M-1,2:M-1)=eye(M-2); %(corner derivatives eliminated)
A(2:M-1,2+M:M+M-1)=-eye(M-2);
topind=[2:M-1:2:M-1];%topindices
%bottom boundary:
A(N*M-M+2:N*M-1,N*M-M+2:N*M-1)=-eye(M-2);
A(N*M-M+2:N*M-1,N*M-M+2-M:N*M-M-1)=eye(M-2);
botind=[N*M-M+2:N*M-1];
%left boundary:
A(M+1:M:(N*M-2*M+1),M+1:M:(N*M-2*M+1))=-eye(N-2);
A(M+1:M:(N*M-2*M+1),M+1+1:M:(N*M-2*M+1+1))=eye(N-2);
leftind=[M+1:M:(N*M-2*M+1)];
%right boundary:
A(2*M:M:(N*M-M),2*M:M:(N*M-M))=eye(N-2);
A(2*M:M:(N*M-M),2*M-1:M:(N*M-M-1))=-eye(N-2);
rightind=[2*M:M:(N*M-M)];
%%Now do corner (diagonal derivatives):
A(1,1)=1;A(1,M+2)=-1;
A(M,M)=1;A(M,2*M-1)=-1;
A(N*M-M+1,N*M-M+1)=-1;A(N*M-M+1,N*M-2*M+2)=1;
A(N*M,N*M)=-1;A(N*M,N*M-M-1)=1;
cornind=[1,M,N*M-M+1,N*M];


%%now set up right hand-side: bc's
% are set and interior source values are set:
b=zeros(M*N,1);
b(2:M-1,1)=zeros(M-2,1);%top bcs DOES NOT INCLUDE CORNERS
b(N*M-M+2:N*M-1,1)=zeros(M-2,1); %bottom bc - see top bc
    %%to satisfy condition of no net "mass" flux into box
b(M+1:M:(N*M-2*M+1),1)=zeros(N-2,1); %left boundary conds.
b(2*M:M:(N*M-M),1)=zeros(N-2,1); %right bcs:
b(1,1)=0;b(M,1)=0; %upper cornersb(N*M-M+1,1)=0;b(N*M,1)=0;%lower corners
%%interior sources:
for j9=2:N-1
  for i9=2:M-1
  k9=(j9-1)*M+i9; %interior grid no. again
  b(k9)= 0*(Deltax)^2; %%note appearance of Deltax. All 0 here.
%%to make Neumann problem consistent - if have pure
      %%flux conditions can't specify interior sources
      %%independently
       end, %end i9 loop
        end, %end j9 loop
b(13)=1*Deltax^2;

%%make up a set of all boundary indices:
bdyind=sort([topind,botind,rightind,leftind,cornind]);
