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


function [P,gtime] =statogram(idata,xdata,deltat,iplot,jump,dt)

idata=idata(:)';
xdata=xdata(:)';

if isempty(xdata), xdata = (1:length(idata))*deltat; end;

count = 0;
for place = 1:jump:length(idata)-dt 
  count  = count + 1;
  pidata = idata(place:place+dt);
  gtime(count) = round(mean(xdata(place:place+dt)));
  P.m(count) = mean(pidata);
  P.v(count) = var(pidata);
  P.s(count) = skewness(pidata);
end;  %Stop program loop

