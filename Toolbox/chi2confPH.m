%function [ci]=chi2confPH(cl,v);
%
%conifidence interval for a chi-squared distribution.
%
%ci confidence interval
%cl confidence level
%v  degrees of freedom


function [ci]=chi2confPH(cl,v);

ci=[chi2inv((1-cl)/2,v)/v chi2inv(1-(1-cl)/2,v)/v];

%ci2=[(1-2./(9*v)-1.96*sqrt(2./(9*v))).^3 (1-2./(9*v)+1.96*sqrt(2./(9*v))).^3];




