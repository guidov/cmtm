%This demeans and normalizes a given input function to unit std.   
%NaNs are ignored and left as are.
%
%[odata]=normPH(idata);

function [odata]=normPH(idata);

%pl=find(isnan(idata)==1); 
%if length(pl)>0;
%display(['There were the following # of NaNs: ',num2str(length(pl))]);
%idata(pl)=nanmean(idata);
%end;
odata=idata-nanmean(idata);
odata=odata/nanstd(odata);