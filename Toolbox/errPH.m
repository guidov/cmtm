%plots error bars.
%
%x= x-axis location
%y= y-axis location
%r= one sided error bar size
%dy=length of one sided y-bar
%m=marker
%s=markersize (last 3 values are optional
%
%function []=errPH(x,y,r,dy,m,s); 


function []=errPH(x,y,r,dy,m,s); 


if nargin<4, dy=.05; end;
if nargin<5, m='k.'; end;
hold on;
for ct=1:length(r);
  h(ct)=plot(x(ct),y(ct),m);
  plot([x(ct)-r(ct) x(ct)+r(ct)],[y(ct) y(ct)],m(1));
  plot([x(ct)-r(ct) x(ct)-r(ct)],[y(ct)-dy y(ct)+dy],m(1));
  plot([x(ct)+r(ct) x(ct)+r(ct)],[y(ct)-dy y(ct)+dy],m(1));
end;

if nargin==6, 
  set(h,'markersize',s);
end;
