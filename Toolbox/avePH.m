%This function averages the input data by N number of points
%The returned data fix(length(Idata)/N) long.
%
% function [Odata] = avePH(Idata,N);


function [Odata] = avePH(Idata,N);

ct2 = 0;
long = length(Idata);
for ct = 1:N:long;
ct2 = ct2+1;
if ct+N<=long; Odata(ct2) = mean(Idata(ct:ct+N-1));
else Odata(ct2)= mean(Idata(ct:long)); end;
end;

