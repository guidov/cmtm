%This makes a new yaxis plot.
%
%x =x-ordinate
%y1=true y-ordinate
%y2=plotted y-ordinate (often de-meaned or the like)
%Nt =number of sig figs to use, or if > 1# tick marks to use
%Ns =number of sig figs to use.
%pos: [] or 1 =left, 2 = right
%
%function []=xaxPH(x,y0,Nt,Ns);

function []=xaxPH(x,y0,dy,Nt,Ns);

dx     = (max(x)-min(x))/200;
yoff   = dy/200;
if length(Nt)==1,
xtick  = linspace(min(x),max(x),Nt);
else,
xtick = Nt;
end;

xtick  = round(xtick);

  plot([min(x) max(x)],[y0 y0],'k');

for ct=1:length(xtick);
  plot([xtick(ct) xtick(ct)],[y0 y0+3*yoff],'k');
  text(xtick(ct)+3*dx,y0-12*yoff,num2str(xtick(ct),Ns));
end;

