%plots error bars.
%
%x= x-axis location
%y= y-axis location
%r= one sided error bar size
%dx=length of one sided y-bar
%m=marker
%s=markersize (last 3 values are optional
%
%function []=errPH(x,y,r,dx,m,s);

function []=erryPH(x,y,r,dx,m,s); 


if nargin<4, dx=.05; end;
if nargin<5, m='k.'; end;
hold on;
for ct=1:length(r);
  h(ct)=plot(x(ct),y(ct),m);
  set(h(ct),'markerfacecolor',m(1));
  if nargin==6, set(h(ct),'markersize',s); end;
  plot([x(ct) x(ct)],[y(ct)-r(ct) y(ct)+r(ct)],m(1));
  h=plot([x(ct)-dx x(ct)+dx],[y(ct)-r(ct) y(ct)-r(ct)],m(1));
  h=plot([x(ct)-dx x(ct)+dx],[y(ct)+r(ct) y(ct)+r(ct)],m(1));
end;

