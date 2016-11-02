%Same as detrend but forces line to pass through the origin.  

function [y0]=detrend0(x,y);

Pf=polyfit0(x,y,1);
Pf=[Pf 0];
Pv=polyval(Pf,x);
y0=y-Pv;
