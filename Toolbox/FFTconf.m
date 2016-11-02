%This checks what the appropriate degrees of freedom are for 
%different scenarios of zero padding and interpolating.
%
%function [lu] = FFTconf(long, pad, level, runs, qplot); 
%
%lu   = upper confidence level
%long = length of the record
%pad  = length of zero padding used
%level= what confidence interval
%runs = number of runs to average


function [lu] = FFTconf(long, pad, level, runs, qplot); 

randn('state',sum(100*clock))
for ct=1:runs;
rr   = randn(1,long); 
rr   = rr-mean(rr);
rr   = rr/std(rr);
[freq P] = fftPH(rr,1,[],pad);
SP       = sort(P);
lu(ct)  = SP(round(level*length(P)));
end;
lu = mean(lu);

if qplot>0;
  figure(qplot); clf; hold on;
  plot(freq,P,'.');
  plot([min(freq) max(freq)],[lu lu],'r--');
  title('Confidence level applied to a sample spectrum');
  ylabel('Power');
  xlabel('Frequency'); 
end;
