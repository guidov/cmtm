%This Calculates the bi-coherence.
%Adapted from the bispectrum routine in Ice Ages and Astronomical Causes
%p. 293, Muller and MacDonald, 2000.
%
% [bispec, freq]=bicohPH(y,dt,qplot,window,fmax)
%
%y      = vector to form bicoherence of.  
%dt     = unit step interval
%qplot  = 1 for plot, 0 for not. 
%window = the default is a hanning window: boxcar and hamming will work.
%fmax   = max freq to calculate up to, 1/17 is the default. 
%dim    = dimension of smoother used in smooth2PH.m (default=3);

function [Bispec, freq]= bicohPH(y,dt,qplot,window,fmax,sm);

if nargin<6, sm=3; end;

A=size(y); if A(1)<A(2), y=y'; end;

if ~exist('fmax'), fmax=.06; end;

pad=length(y);
y=y-mean(y);

%Window the data
if length(window)==0, window='hanning'; end;
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
B  = zeros(nfreq);

%Calculate with a normalization.
for i = 1:nfreq,
for j = 1:nfreq,
    B(i,j)=fy(i)*fy(j)*conj(fy(i+j-1))./(abs(fy(i))*abs(fy(j))*abs(fy(i+j-1)));
end; end;

%Smooth the bispectrum
%filt = [.0 .5 .0; .5 1 .5; .0 .5 .0];
BC    = abs(smooth2PH(B,sm));
%NC    = smooth2PH(abs(B),sm);
Bispec = BC(3:nfreq,3:nfreq); 
%Bispec = BC(3:nfreq,3:nfreq)./NC(3:nfreq,3:nfreq);
freq   = freq(3:nfreq);

if qplot==1; 
%figure(1); clf;
hold on;
Bispec(end,end)=.6;
Bispec(end,end-1)=1;
contourf(freq,freq,Bispec,[.4:.05:1]);
axis tight;
%xlabel('Frequency');
%ylabel('Frequency');
%title('Bispectrum');
colorbar('vert');
end;


