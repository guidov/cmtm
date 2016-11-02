%This funciton finds the Q value of a peak - for use with Lcore
%Difference from qPHa is that this only marks stars at the peaks
%Difference from qPHb is this gets rid of freqs with Milank. tendencies
%
%Format of the function is this:
%
%  [xdata, Pfreq] = qPHb(ydata,xdata,block,freq)
%
% xdata is returned with Milank. freqs set to NaN
% Pband is the sum of the energy at the Milankovitch Bands
% block = data points to average over in looking for peak edges.
% freq  = specify freq to automatically assing Q values to


function [Nxdata,Pband] = qPHb(ydata,xdata,block,freq)

if isempty(xdata), 
  df = 1/(length(ydata));
  xdata = df:df:df*length(P)
end;
Nxdata = xdata;
for ct = 1:length(freq);
   x = freq(ct);
   [i, px] = min(abs(xdata-x));
if (px>block) & (px<length(ydata)-block); 
   [i, p] = max(ydata(px-block:px+block));
   p = p+px-(block+1);
else;
p=px;
if p==length(ydata); p=px-1; end;
if p==1;             p=2;    end; 
end;
   while ydata(p)<ydata(p+1) & p+2<length(ydata); p = p+1; end;
   while ydata(p)<ydata(p-1) & p-1>1;             p = p-1; end;
   if p-block < 1, blockl = p-1; else blockl = block; end; 
   if p+block > length(ydata), blockr = length(ydata)-p; else blockr = block; end;

   [i, pl] = min(ydata(p-blockl:p));
   pl = pl+p-blockl-1;
   while ydata(pl)<ydata(pl+1) & pl>1; pl=pl-1; end;
   [i, pr] = min(ydata(p:p+blockr));
   pr = pr+p-1;
if pr+1>length(ydata); pr=pr-1; end;
   while ydata(pr+1)<ydata(pr) & pr<length(ydata)-1; pr=pr+1; end;

   mhi = (2*ydata(p)+ydata(pl)+ydata(pr))/4;
   [i, lm] = min(abs(ydata(pl:p-1)-mhi));
   [i, rm] = min(abs(ydata(p+1:pr)-mhi));
   lm = lm + pl-1;
   rm = rm + p;
      
   Nxdata(pl:pr) = NaN;
   Pband(ct) = sum(ydata(pl:pr));   

   %some plotting of the results
   hold on
   plot(xdata(p),ydata(p),'r*');
   H = text(xdata(p),ydata(p),[' ',num2str(1/xdata(p),3)]); set(H,'fontsize',10); 

   
   clear pry ply i pr pl i lm rm p hi mhi width;
end;  %stop freq for statement.




