%This function computes the cross correlation (r^2) at zero lag for
%two input records.  Records should be sampled at the same 
%unifrom intervals.
%
%records y1 and y2, default is r^2, but if type==1 then returns r.
%
%function [XC]=xcPH(y1,y2,type,demean);
%
%y1:     1st record 
%y2:     2nd record
%type:   0=squared cross-correlation, 1=cross-correlation       (default=0)
%demean: 0=leave the means of y1 and y2 alone, 1=take mean out. (default=1)
%window  0=no, 1=yes, use a Hanning window (default=0);

function [XC]=xcPH(y1,y2,type,demean,window);

y1=y1(:); y2=y2(:);

if nargin<3, type=0; end;
if nargin<4, demean=1; end;
if nargin<5, window=0; end;

pl=find(isnan(y1)==0 & isnan(y2)==0);
y1=y1(pl); y2=y2(pl);

if demean==1,
y1=y1-mean(y1);
y2=y2-mean(y2);
end;

if window==1,
  y1=y1.*hanning(length(y1));
  y2=y2.*hanning(length(y2));
end;

if isreal(y1) & isreal(y2),
XC= sum(y1.*y2)/sqrt(sum(y1.*y1)*sum(y2.*y2));
else, 
XC= sum(y1.*conj(y2))/sqrt(sum(y1.*conj(y1))*sum(y2.*conj(y2)));
end;

if type~=1, XC=XC^2; end; 
