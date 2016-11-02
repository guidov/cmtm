% Creates a spectrogram of a time series by taking an interval of the time
% series with length 'dt' and sliding it foward by interval 'jump' along
% the length of the time series.
%
% Input should read as follows:
%
% [P,freq,gtime] =specPH(idata,xdata,deltat,iplot,jump,dt)
%
%
% idata  = input for spectrogram
% xdata  = data put on x axis, i.e. mean time associated with segment
% deltat = time between samples
% iplot  = do you want to plot the results
% jump   = move chunk over this many spaces in the data
% dt = length of segment to take power spectrum of
%
% The power spectrum of each chunk is computed using pmtmPH.  The pmtmPH
% function automatically uses 4 window and all input data Function display
% the mean, variance, and variance times deltat for each segment spectrally
% analyzed.  Further it will get rid of NaNs and INFs in your data and tell
% you how many it got rid of.  This is a computationally intensive function


function [P,freq,gtime] =spectrogram(idata,xdata,deltat,iplot,jump,dt,minf,maxf)

idata=idata(:)';
xdata=xdata(:)';

if nargin<8, maxf=1/deltat; end;
if nargin<7, minf=0; end;

if isempty(xdata), xdata = (1:length(idata))*deltat; end;


idata=normPH(idata);
count = 0;
for place = 1:jump:length(idata)-dt 
  count  = count + 1;
  pidata = idata(place:place+dt);
  %pidata = hanning(length(pidata))'.*pidata;
  %pidata = tukeywin(length(pidata),.5)'.*pidata;
  gtime(count) = round(mean(xdata(place:place+dt)));
  [Pxx,freq,Pxxc] = pmtmPH(detrend(pidata),deltat,0,0);
  %[freq,Pxx] = fftPH(pidata,deltat,[]);
  Pxx=Pxx';
  %Pxx=Pxx.*freq;
  P(1:length(Pxx),count) = Pxx;
end;  %Stop program loop

gtime(length(gtime)+1) = gtime(length(gtime));
P(:,length(P(1,:))+1) = P(:,length(P(1,:)));

P=log10(P/max(max(P)));

%P=smooth2PH(P,3); 
min(min(P)),
max(max(P)),

if (iplot) == 1
  pl=find(P<-2); P(pl)=-2;
  %pl=find(P<-3.1); P(pl)=-3.1;
  %pl=find(P>-0.4); P(pl)=-0.4;
  %P=P*(max(max(P))-min(min(P)));
  pl = find(freq>minf & freq<maxf);
  pcolorPH(gtime,freq(pl),P(pl,:));
  shading flat;  colorbar('vert');
  axis tight;   
  h=ylabel('frequency'); font(h,14);  font(gca,12);
  set(gca,'ydir','reverse');
end;
