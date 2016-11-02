%Makes an arrow
%x,y location
%dir 1 (right), 2(left), 3 (up), or 4 (down)
%color
%
%function arrowPH(x,y,dir,col,size); 

function arrowPH(x,y,dir,col,size); 

hold on;

if nargin==0, x=0, y=0; end;
if nargin<3, dir=1; end;
if nargin<4, col='k'; end;
if nargin<5, size=1; end;

if dir==1 | dir==2,
  xx = [-.25 -.25 .25];
  yy = [.15 -.15 0];
  if dir==2, xx=-xx; end;
end;

if dir==3| dir==4,
  yy = [-.25 -.25 .25];
  xx = [.15 -.15 0];
  if dir==4, yy=-yy; end;
end;

h=fill(size*xx+x,size*yy+y,'k');  
set(h,'edgecolor',col); 
set(h,'facecolor',col); 

%h=plot([0 1],[0 0],'k');
%set(h,'color',col);
%set(h,'linewidth',2); 


