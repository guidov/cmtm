%Calculates the bispectrum.
%Taken from Ice Ages and Astronomical Causes p. 293, Muller and MacDonald, 2000.
%
% bispectrum(y,dt,qplot,window,fmax)
%
%The default is a hanning window: boxcar and hamming will work.
%fmax = max freq to calculate up to. 

function [Bispec, freq]= bispectrum(y,dt,qplot,window,fmax);

if ~exist('fmax'), fmax=.06; end;

pad=2^13;
y=y-mean(y);

%Window the data
if ~exist('window'), window='hanning'; end;
eval(['y=y.*',window,'(length(y));']);
fy = fft(y,pad);

Nyquist = .5/dt;
df    =   1/(pad*dt);
freq  = 0:df:Nyquist;

%Discard low freqs. due to padding.
pl = find(freq>1/(length(y)*dt));
freq = freq(pl); fy=fy(pl);

%Find number of freqs to look at. 
nfreq=fix(fmax/df);
B     = zeros(nfreq);

for i = 1:nfreq,
for j = 1:nfreq,
    B(i,j)=fy(i).*fy(j).*conj(fy(i+j-1));
end; end;

%Smooth the bispectrum
filt = [.0 .5 .0; .5 1 .5; .0 .5 .0];
BC   = abs(conv2(filt,B)).^2;
Bispec = BC(2:nfreq,2:nfreq);
freq   = freq(1:nfreq-1);

Bispec=Bispec/mean(mean(Bispec));

if qplot==1; 
figure(1); clf;
contour(freq,freq,log10(Bispec),10);
axis tight;
xlabel('Frequency');
ylabel('Frequency');
title('Bispectrum');
colorbar('vert');
end;

