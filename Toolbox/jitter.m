%This function makes a depth profile using SARs with a spectrum as derived in tsar.m 
%Must have 2 control points only.
%Uses 1/(1/100+f) frequency scaling.
%
%function [Dep,Lt]=MakeDepSimple(chatter,Tt);
%
%Dep     = Depth profile
%Lt      = Linear Time associated with the depth profile and NCP
%jitter  = variance of the SAR divided by the mean SAR^2.
%Tt      = true time
%fo      = decorelation frequency, default is 1/100; 

function [Dep,Lt]=jitter(jitter,t,fo);

if nargin<3, fo=1/100; end;
dt  = mean(diff(t));
df  = 1/(max(t)-min(t));
f   = df:df:1/(2*dt);
Filter=sqrt(1./(fo^2+f.^2));
%Filter = sqrt(f.^(-1));  %number is the power law of S.
Filter = [max(Filter) Filter fliplr(Filter)];
x   =randn(length(t),1)'; x=normPH(x);
fx  =fft(x);
ffx =fx.*Filter;
SAR  =normPH(real(ifft(ffx)));

SAR=SAR*sqrt(jitter)+1;
pl=find(SAR<=0);
SAR(pl)=.0001;

Dep=cumsum(SAR);
Dep=Dep-min(Dep);
Lt=interpPH([Dep(1) Dep(end)],[t(1) t(end)],Dep);           





