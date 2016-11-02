%This checks what the appropriate degrees of freedom are for different scenarios of zero padding and interpolating.
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
long = 2^11;
pad  = 0;
ct = 0;
randn('state',sum(100*clock))
rr   = randn(1,long); 
rr   = rr-mean(rr);
rr   = rr/std(rr);
[freq P] = fftPH(rr,1,[],pad);
SP       = sort(P);
S95     = SP(round(.95*length(SP)));
S025     = SP(round(.025*length(SP)));  
lsp      = linspace(1/length(P),1,length(P));

figure(1); clf; hold on;
plot(lsp,SP);
plot([min(lsp) max(lsp)],[S95 S95],'k:');
plot([.95 .95],[min(P) max(P)],'k:');
