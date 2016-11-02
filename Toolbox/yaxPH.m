%This makes a new yaxis plot.
%
%x =x-ordinate
%y1=true y-ordinate
%y2=plotted y-ordinate (often de-meaned or the like)
%Nt =number of tick marks 
%Ns =number of sig figs to use
%pos: [] or 1 =left, 2 = right
%
%function [h]=yaxPH(x,y1,y2,Nt,Ns,pos);

function [h]=yaxPH(x,y1,y2,Nt,Ns,pos,hs);


dx     = max(x)-min(x);
xo     = dx/200;
ytick  = linspace(min(y2),max(y2),Nt);
ylabel = linspace(min(y1),max(y1),Nt);

if nargin<6 | pos==1,
  xc=min(x)-dx/50;
  plot([xc xc],[min(y2) max(y2)],'k');
  for ct=1:Nt;
    plot([xc-.75 xc],[ytick(ct) ytick(ct)],'k');
    h=text(xc+3,ytick(ct),num2str(ylabel(ct),Ns));
    
    if nargin>6,
      font(h,hs); 
    end;

  end;
else,
  xc=max(x)+dx/100;
  plot([xc xc],[min(y2) max(y2)],'k');
  for ct=1:Nt;
    plot([xc xc+.75],[ytick(ct) ytick(ct)],'k');
    h=text(xc-1.1,ytick(ct),num2str(ylabel(ct),Ns));
  
    if nargin>6,
      font(h,hs); 
    end;
  end;
end;
  
