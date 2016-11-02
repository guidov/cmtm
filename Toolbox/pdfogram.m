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


function [P,gtime] =pdfogram(idata,xdata,deltat,iplot,jump,dt,bins)

idata=idata(:)';
xdata=xdata(:)';

if nargin<7, bins=10; end;

if isempty(xdata), xdata = (1:length(idata))*deltat; end;

count = 0;
edges=linspace(mean(idata)-3*std(idata),mean(idata)+3*std(idata),bins);
for place = 1:jump:length(idata)-dt+1 
  count  = count + 1;
  pidata = idata(place:place+dt);
  gtime(count) = round(mean(xdata(place:place+dt)));
  P(:,count) = histc(pidata,edges);
  P(:,count)=P(:,count)/sum(P(:,count));
end;  %Stop program loop

gtime(length(gtime)+1) = gtime(length(gtime));
P(:,length(P(1,:))+1) = P(:,length(P(1,:)));

if (iplot) == 1
  pcolorPH(gtime,edges,P);
  shading flat;  colorbar('vert');
  axis tight;   
  h=ylabel('values'); font(h,14);  font(gca,12);
end;