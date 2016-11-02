%Interpolating in a rational way.
%
%qext=1 for extrapolation of mean rate 
%qest=0 for NaN outside of x1 limits.
%
%function [y2]=interpPH(x1,y1,x2,qext);

function [y2]=interpPH(x1,y1,x2,qext);

if length(x1)~=length(y1), 
  disp('input vectors are not the same length'); 
  return;
end;

if nargin<3,  disp('Not enough input'); return; end;
if nargin==3, qext=0; end;


[x1 j]=sort(x1); y1=y1(j); 
[a b]=size(x1); if a>b, x1=x1'; end;
[a b]=size(y1); if a>b, y1=y1'; end;


for ct=1:length(x2);
     pl1=find(x1>=x2(ct));
     if length(pl1)==0, 
       y2(ct)=NaN; 
     else,
       pl=pl1(1);
       if x2(ct)==x1(pl), 
	 y2(ct)=y1(pl); 
       else,
         if pl==1, 
	   y2(ct)=NaN;      
	 else,
           y2(ct)=(x2(ct)-x1(pl-1))*(y1(pl)-y1(pl-1))/(x1(pl)-x1(pl-1))+y1(pl-1);
         end;
       end;
     end;
end;

if qext==1,
  pl2=find(isnan(y2)==1);
  if length(pl2)>0,
    pl=find(isnan(y1)==0 & isnan(x1)==0);
    S=(y1(pl(1))-y1(pl(end)))/(x1(pl(1))-x1(pl(end)));
    b=y1(pl(1))-S*x1(pl(1));
    y2(pl2)=x2(pl2)*S+b;
  end;
end;

