%Makes noise with a given power-law.
%
%function [noise]=rednoise(t,plaw,fo);
%
%t    = time.
%plaw = power law spectra follows for freqs above fo.
%fo   = decorelation frequency, default is 1/infinity.

function [noise]=rednoise(t,plaw,fo);

if nargin<2, plaw=2;   end;
if nargin<3, fo=0; end;

dt     = mean(diff(t));
df     = 1/(max(t)-min(t));
f      = df:df:1/(2*dt);
Filter = sqrt(1./(fo^plaw+f.^plaw));
Filter = [max(Filter) Filter fliplr(Filter)];
Filter(Filter==0)=1; 
x      = randn(length(t),1)'; x=normPH(x);
fx     = fft(x);
ffx    = fx.*Filter;
noise  = normPH(real(ifft(ffx)));





