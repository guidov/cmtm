%Confidence interval using a monte carlos method
%for use with a padded periodogram.
%A sorted average power vector is generated to compare ones values against.   
%
%function [Confid]=ChiPH(M,Pad);
%
%M   = number of measurements;
%Pad = How much zero padding as the exponent of two: always 13 for now.

%function [Confid]=ChiPH2(Pval,M,Pad);
clear;
M = 2^8;
Pad = 13;
Pval = 2;
load Confid_Table;

%Add to the table if type does not yet exist.
pl2=find(Type(:,1)==M);
if isempty(pl2); pl2=NaN; end;
pl3=find(Type(:,2)==Pad);
if isempty(pl2); pl3=0;   end;
[i pl] = find(pl2==pl3);
if isempty(pl);
  Type = [Type; M Pad]
  randn('state',sum(100*clock))
  for ct = 1:100;
    [freq P] = fftPH(randn(M,1),1,[],Pad);
    SP(ct,:) = sort(P);
  end;
  MSP = mean(SP);  %Mean sorted power;
Table(:,length(Table(1,:))+1) = MSP';
save Confid_Table Table Type;
else;
MSP = Table(:,pl);
end;

[i pl] = min(abs(MSP-Pval));
Conf   = pl/length(MSP(:,1));

