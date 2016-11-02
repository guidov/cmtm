%function [Covariance] = covPH(x,y,lags);
%
%The difference between this function and xcov, is that it normalizes each
%lag by the sqrt(sum((x-xm)^2)*sum((y-ym)^2)).  This is Pearsons r, r^2 gives the
%degree by which the two records covary (how much of the variance in y is accounted
%for the by the variance in x.  

function [xc] = covPH(xi,yi,lags);

Ni = length(xi);
for mv = 1:(2*lags+1);
  lag = mv-(lags+1);
  if lag<=0;
    x = xi(1:Ni-abs(lag));
    y = yi(abs(lag)+1:Ni);
  else
    x = xi(lag+1:Ni);
    y = yi(1:Ni-lag);
  end;
xc(mv)  = xcov(x,y,0)/sqrt(sum((x-mean(x)).^2)*sum((y-mean(y)).^2));
end;
