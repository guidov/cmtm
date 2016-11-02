%Chi squared confidence level table look-up function.
%Using the tabel made with ChiTable.m, this function looks up 
%the confidence level for a given normalized value and its 
%associated degrees of freedom assuming a chi-squared distribution.  
%
%function [Confid]=ChiPH(NL,df);
%
%NL = Normalized level = Peak Value divided by local background noise level.
%df = integer degrees of freedom.
%Confid = the confidence you can have in this particular peak.

function [Confid]=ChiPH(NL,df);

load ChiTable;

level  = level(:,df);           %Only use column with given degrees of freedom.
Confid = interp1(level,P,NL);  %Interpolate confidence to given normalized level.
if Confid>.99; Confid=1; end;
if isnan(Confid); Confid=0; end;
if NL<=1; Confid=0; end;


