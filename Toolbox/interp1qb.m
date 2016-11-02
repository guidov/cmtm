function yi=interp1qb(x,y,xi)
%   INTERP1Q Quick 1-D linear interpolation..
%
%   F=INTERP1Q(X,Y,XI) returns the value of the 1-D function Y at the
%   points XI using linear interpolation. Length(F)=length(XI).
%   The vector X specifies the coordinates of the underlying interval.
%
%   NaN's are returned for values of XI outside the coordinates in X.
%   X must be a monotonically increasing.


    x=x(:); y=y(:); xi=xi(:);
    pl=find(isnan(x)==0 & isnan(y)==0);
    x=x(pl); y=y(pl);
    
   [xxi,k] = sort(xi);
   [dum,j] = sort([x; xxi]);
   r(j)    = 1:length(j);
   r       = r(length(x)+1:end)-(1:length(xxi));
   r(k)    = r;
   r(xi==x(end))=length(x)-1;
   pl      = find((r>0) & (r<length(x)));
   rpl     = r(pl);
   u       = (xi(pl)-x(rpl))./(x(rpl+1)-x(rpl));
   yi      = repmat(NaN,length(xxi),1);
   yi(pl)  = y(rpl,:)+(y(rpl+1,:)-y(rpl,:)).*u;
