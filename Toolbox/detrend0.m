%Same as detrend but forces line to pass through the origin.  
%if x is included, then non-evenly space data may be detreneded.
%if not inlcuded then its assumed the y begins at 0 and is evenly spaced.
%
%function [y0]=detrend0(y,x);

function [y0]=detrend0(y,x);

if ~exist('x'); x=0:length(y)-1; end;

if length(y(:,1))~=length(x(:,1)); y=y'; end;

Pf=polyfit0(x,y,1);
Pf=[Pf 0];
Pv=polyval(Pf,x);
y0=y-Pv;
