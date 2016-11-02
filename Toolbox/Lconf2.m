%This checks what the appropriate degrees of freedom are for 
%different scenarios of zero padding and interpolating.
%
%function [lu ll] = FFTconf(long, pad, level, runs); 
%
%lu   = upper confidence level
%ll   = lower level
%long = length of the record
%pad  = length of zero padding used
%level= what confidence interval
%runs = number of runs to average

clear;
long = 2^14;
pad  = 0;
ct = 0;
nw = 3;
for a= .05;
  ct = ct+1;
randn('state',sum(100*clock))
rr   = randn(1,long); 
rr   = rr-mean(rr);
rr   = rr/std(rr);
[P Pc freq] = pmtm(rr,nw);
SP       = sort(P);
S95      = SP(round((1-a)*length(SP)));
e = chi2conf((1-a),2*nw-1);
%e = e*mean(P);
level(ct) = a;
ratio(ct) = e(1)/S95;
lrdiff(ct)    = e(1)-S95;
end;

clf; hold on;
loglog(level,lrdiff);
