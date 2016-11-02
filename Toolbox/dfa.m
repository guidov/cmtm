%An implementation of Detrended Fluctuation Analysis
%
%See Peng et al, "Mosaic organization of DNA nucleotides"
%Physical Review E, v49 no2 Feb 1994.
%
%   dt: sampling interval
%   y:  sample values
%   qplot: do you want a plot of the output? (1=yes, 0 or omit = No)
%
%   F:  average detrended RMS variability as a function of scale, l. 
%   L:  interval length
%   slope:  slope of log10(F)/log10(L)
%
%function [F,L,slope]=dfa(y,l,qplot);

function [F,L,slope]=dfa(y,dt,qplot);

dt=abs(dt);
y=cumsum(y);

y=y(:);
N=length(y);
L=4:N/4; 

if nargin<2,
  dt=1;
end;
if nargin<3,
  qplot=0;
end;

for ct=1:length(L),
    %divide into n columns with L rows
    n=floor(N/L(ct));
    M=reshape(y(1:n*L(ct)),L(ct),n);

    %detrend each column
    M=detrend(M); 
    
    %rms of the residuals
    F(ct)= ( sqrt(sum(sum(M.^2))/N) );
end;
L=L*dt;
pl=find(L>0 & F>0);
[P,S]=polyfit(log10(L(pl)),log10(F(pl)),1);
slope=P(1);
if qplot==1,
  %L=log10(L); F=log10(F);
  %figure(1); clf; hold on;
  plot(L,F,'ro'); 
  plot(L,10.^polyval(P,log10(L)),'k'); 
  text(L(1),F(end),num2str(P(1)));
  axis tight; logPH;
end;