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

function [Bispec, freq]= bicohxcPH(x,y,dt,qplot,window,fmax);

A=size(y); if A(1)<A(2), y=y'; end;
A=size(x); if A(1)<A(2), x=x'; end;

if nargin<4, qplot=0; end;
if nargin<5, window='hanning'; end;
if nargin<6, fmax=.06; end;


pad=length(y);
%y=y-mean(y);
%x=x-mean(x);

%Window the data

eval(['y=y.*',window,'(length(y));']);
fy = fft(y,pad);
eval(['x=x.*',window,'(length(x));']);
fx = fft(x,pad);


Nyquist = .5/dt;
df      =  1/(pad*dt);
freq    =  0:df:Nyquist;

%Discard low freqs. due to padding.
pl   = find(freq>1/(length(y)*dt));
freq = freq(pl); 
fy=fy(pl);
fx=fx(pl);

%Find number of freqs to look at. 
nfreq=fix(fmax/df);
B=zeros(nfreq);

%Calculate with a normalization.
for i = 1:nfreq,
for j = 1:nfreq,
    B(i,j)=fx(i)*fx(j)*conj(fy(i+j-1)); %./(abs(fy(i))*abs(fy(j))*abs(fy(i+j-1)));
end; end;

%Smooth the bispectrum
big=13;
cut=(big-1)/2;
A=ones(big,big);
BC = abs(conv2(A,B));
BC = BC(cut:end-cut,cut:end-cut);

%BC    = abs(smooth2PH(B));
NC = conv2(A,abs(B));
NC = NC(cut:end-cut,cut:end-cut);
Bispec = BC(1:nfreq,1:nfreq)./NC(1:nfreq,1:nfreq);
freq   = freq(1:nfreq);

Bispec = Bispec(cut:end-cut,cut:end-cut);
freq   = freq(cut:end-cut);



if qplot==1; 
pl=find(Bispec<.3); 
Bispec(pl)=NaN;
pcolorPH(freq,freq,Bispec); shading flat;
axis tight;
xlabel('Frequency (1/Kyr)');
ylabel('Frequency (1/Kyr)');
colorbar('horiz');
end;
