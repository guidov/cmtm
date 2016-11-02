% function [Plaw,Pval,freq] = plawPH(P,freq,Fmin,Fmax,qplot);
%
%  defaults:
%  qplot  = 0 (no plot)
%  Fmax   = nyquist freq
%  Fmin   = 0
%  freq   = based on dt=1 (aproximate)
%
%see also statPH for another implementation of this

function [Plaw,Pval,freq] = plawPH(P,freq,Fmin,Fmax,qplot);


if nargin<5;  qplot=0; end;
if ~exist('Fmax') | length(Fmax)==0, max(freq)=0; end;
if ~exist('Fmin') | length(Fmin)==0, min(freq)=0; end;
if ~exist('freq') | length(freq)==0, freq=(1:length(P)/2)/length(P); end;
if nargin<1;  help plawPH, return, end;
if (size(P)~=size(freq)) & (size(P)==size(freq'));
  freq=freq';
end

pl   = find( (freq > Fmin) & (freq<=Fmax) ); freq=freq(pl); P=P(pl);
Plaw = polyfit(log10(freq),log10(P),1 );
Pval = polyval(Plaw,log10(freq)); Pval=10.^Pval;
Plaw = Plaw(1);

if qplot,
  hold on;
  h=plot(freq,Pval,'k'); set(h,'linewidth',1.2);
  h=text(mean(freq)/1.75,mean(Pval),['      S=' num2str(round(10*Plaw)/10)]); font(h,14);
end;


