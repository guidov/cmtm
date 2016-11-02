%This function uses Spearmans Rank-Order Correlation Coefficient 
%Ref. p488 Numerical Recipes, Press et al.  
%By ranking the input values the we now know the probability distribution
%and can better estimate the significance of the correlation factor.
%Rs is distibuted approximately as a Students distribution with N-2 
%degrees of freedom.
%
%function [Rs Sig] = Ncorr(x,y,lags);
%
%Rs is the correlation coefficient and Sig is the significance level of that correlation
%based a scaling of Rs that has a distribution close to a students T distritbution.

function [Rs, Sig] = Ncorr(x,y,lags);

%Call the subfunction Crank to sort the data
xr = crank(x);
yr = crank(y);

xm = mean(xr);
ym = mean(yr);

%for ct=-lags:lags;
Rs = sum((xr-xm).*(yr-ym))/(sqrt(sum((xr-xm).^2))*sqrt(sum((yr-ym).^2)));
t  = Rs*sqrt((length(xr)-2)/(1-Rs^2));

%The probability the null hypothesis that x and y are uncorrelated is false
%this is the significance level calculated with the cumulative distribution function
%of a students T distribution suing the betainc.m function.  
df  = length(x)-2;  %Number of degrees of freedom

Sig = tp(abs(t),df)-tp(-abs(t),df);

keyboard


%Rank the input data.
function [R] = crank(input);
[S I] = sort(input);

j = 1;
while j < length(S);
if S(j)~=S(j+1); R(j)=j; 
else; 
jstart = j;
while S(j)==S(j+1) & j<length(S); j=j+1; end;
jend   = j;
R(jstart:jend)=(jstart+jend)/2;
end;  %stop the if loop
j = j+1;
end;  %stop the while loop;
if S(j)~=S(j-1); R(j)=j; else R(j)=(jstart+jstop)/2; end;  %Ranking the last data point

keyboard





