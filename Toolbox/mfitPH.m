%This function seeks to fit some sine and cose waves to some data; it should be
%quicker and simpler than fft since it only looks for a few freqs.
%
%function[E] = mfit(y,step,period)
%
%E is the fourier coefficients squared divided by length squared of each freq
%y is the time series to be analyzed
%step is the delta t of y
%period are the ones you are interested at looking at
%
%No trends or means are removed, that is up to you: see detrend and fft

function[E] = mfit(y,step,period)

for scount = 1:length(period)
   a = 0;
   T = period(scount);
   for ycount = 1:length(y)
     a = a + y(ycount)*sin(2*pi*step*ycount/T) + i*y(ycount)*cos(2*pi*step*ycount/T);
  end;
   E(scount) = a.*conj(a)/(length(y)^2);
end;

