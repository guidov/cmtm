%This function uses Spearmans Rank-Order Correlation Coefficient 
%Ref. p488 Numerical Recipes, Press et al.  
%By ranking the input values the we now know the probability distribution
%and can better estimate the significance of the correlation factor.
%Rs is distibuted approximately as a Students distribution with N-2 
%degrees of freedom.
%
%function [xc,Rs,SR,SD] = Ncorr(x,y,lags);
%
%xc = is the normalizde cross correlation / covariance;
%Rs = Correlation of the nonparametric rankings
%SR = The significance of Rs, based on a student T test
%SD = Significance of D based on a normal distribution, and its expected value and variance.


%Hypothesis is that the data sets are uncorellated, this is called a null-hypothesis.
%If we can reject the null hypothesis then we have determined a correlation.
%I report the significance levels as the significane of there being a correlation 
%betweent the two records. 

function [xc, Rs, SR, SD] = Ncorr(xi,yi,lags);

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

%Call the subfunction Crank to sort the data
[xr xf] = crank(x);   %The ranking a number of ties in each group
[yr yf] = crank(y);

xm = mean(xr);
ym = mean(yr);
N  = length(x);

Rs(mv)  = sum((xr-xm).*(yr-ym))/(sqrt(sum((xr-xm).^2))*sqrt(sum((yr-ym).^2)));
t   = Rs(mv)*sqrt((N-2)/(1-Rs(mv)^2));
df  = N-2;  %Number of degrees of freedom
SR(mv)  = tp(abs(t),df)-tp(-abs(t),df);
%The probability the null hypothesis that x and y are uncorrelated is false
%this is the significance level calculated with the cumulative distribution function
%of a students T distribution suing the betainc.m function.  


%Now calculating the sum squared difference of ranks as well; independent check of significance
N   = length(x);
D   = sum((xr-yr).^2);
Dave= 1/6*(N^3-N)-1/12*sum(xf.^3-sf)-1/12*sum(yf.^3-yf);  %Expectation value for null hypothesis
Dvar= 1/36*(N-1)*N^2*(N+1)^2*(1-sum(xf.^3-xf)/(N^3-N))*(1-sum(yf.^3-yf)/(N^3-N));
ZD  = (D-Dave)/sqrt(2*Dvar);
SD(mv)  = erf(abs(ZD));



%And finally the old stand by, cross correlation, normalized as such.  
xc(mv)  = xcov(x,y,0)/sqrt(sum((x-mean(x)).^2)*sum((y-mean(y)).^2));

end;

%keyboard;

%------------------------------------------------------------%
%Rank the input data.
function [Rank, f] = crank(input);
[S I] = sort(input);
f = [];
k = 0;
j = 1;
while j < length(S);
if S(j)~=S(j+1); R(j)=j; 
  else; 
  k = k+1;
  jstart = j;
  while S(j)==S(j+1) & j<length(S)-1; j=j+1; end;
  jend   = j; 
  R(jstart:jend)=(jstart+jend)/2;
  f(k) = jend-jstart+1;    %The number of the ties in this group
end;  %stop the if loop;
j = j+1;
end;  %stop the while loop;
if S(j)~=S(j-1); R(j)=j; else R(j)=(jstart+jend)/2; end;  %Ranking the last data point
Rank(I) = R;    %Putting ranked measurements back in their original order.


