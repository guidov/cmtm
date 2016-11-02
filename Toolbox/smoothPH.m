function [Odata] = smoothPH(Idata,N,filter);

%This function smooths the input data by N number of points
%
%  
% function [Odata] = smoothPH(Idata,N,filter);
%
% N      = number of points to smooth over (rounds up to an odd number)
% filter = shape to convolve with: default is hamming. 
%          Boxcar, hanning, blackman are available.


% make N odd
if rem(N,2)==0, N=N+1; end;

% set defaults
if nargin<3;
  F=hanning(N);
else,
  eval(['F=',filter,'(N);']);
end
    
% normalize area of F to unit area
F=F/sum(F);
Odata=conv(Idata,F)./conv(ones(size(Idata)),F);
Odata(1:floor(N/2))=[];
Odata(end-floor(N/2)+1:end)=[];






