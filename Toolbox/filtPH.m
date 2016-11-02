%This will return a frequency band in the time domain.
%All data will be demeaned.  Method is take fft, set unwanted freqs to zero,
%Take ifft, and return the data to you.
%
%function [odata] = filtPH(idata,deltat,lf,hf,window);
%
%  idata = input data
%  odata = output data
% deltat = time between measurements, neccessary in defining lf and hf appropriately
%  lf/hf = low/high frequency cut off
% window = what to window data with.  Boxcar is default, hanning or hamming works.

function [odata, Pratio] = filtPH(idata,deltat,lf,hf,window); 

A=size(idata); if A(1)<A(2); idata=idata'; end;
idata  = idata-mean(idata);

if nargin<5 | length(window)==0, window='boxcar'; end;
eval(['idata=idata.*',window,'(length(idata));']);


%Take fourier transform
fy     = fft(idata);
Pidata = sum(fy.*conj(fy));

%Make frequency axis, accomodating for the way Matlab outputs fft
long = fix(length(fy)/2);
df   = 1/(length(idata)*deltat);
if length(idata)/2==long,          %Different freq axis if length is even or odd
freq = [0 df:df:df*long (df*long-df):-df:df];
else   freq = [0 df:df:df*long df*long:-df:df]; 
end;


%Find data not in the band, take ifft of the rest
place = [find(freq<=lf) find(freq>=hf)];
fy(place) = 0;

%Applying a window to reduce Gibbs phenomena
%pl = find(freq>=lf & freq<=hf);
%[a split] = max(freq);
%pl1= pl(find(pl<=split));
%pl2= pl(find(pl>split));
%H1 = hamming(length(pl1));
%H2 = hamming(length(pl2));
%fy2 = fy;
%fy2(pl1)=fy(pl1).*H1';
%fy2(pl2)=fy(pl2).*H2';

odata = real(ifft(fy));

%Podata = sum(fy.*conj(fy));
%Pratio = sqrt(Podata/Pidata);
%odata = normPH(odata)*Pratio;
