%This funciton finds the Q value of a peak - designed for use with pmtmPH
%Difference from qPHa is that this only marks stars at the peaks
%
%Format of the function is this:
%
%  [Q, mfreq, Pfreq] = qPHb(ydata,xdata,block,qplot,freq,symbol)
%
% Q     = quality factor
% mfreq = frequency where max power located
% Pfreq = Power contained within peak
% block = data points to average over in looking for peak edges.
% qplot = do you want stuff displayed on the graph; 1 = yes.
% freq  = specify freq to automatically assing Q values to


function [Q,mfreq,Pfreq] = qPHb(ydata,xdata,block,qplot,freq,symbol)


format compact
hold on
if isempty(xdata), 
  df = 1/(length(ydata));
  xdata = df:df:df*length(P)
end;

%mx = mean(abs(xdata));
%my = mean(abs(ydata));


for ct = 1:length(freq);

   x = freq(ct);
   [i, px] = min(abs(xdata-x));
   [i, p] = max(ydata(px-block:px+block));
   p = p+px-(block+1);
   while ydata(p) < ydata(p+1), p = p+1; end;
   while ydata(p) < ydata(p-1), p = p-1; end;
   if p-block < 1, blockl = p-1; else blockl = block; end; 
   if p+block > length(ydata), blockr = length(ydata)-p; else blockr = block; end;

   [i, pl] = min(ydata(p-blockl:p));
   pl = pl+p-blockl-1;
   while (ydata(pl) < ydata(pl+1)) & (pl > 1), pl = pl-1; end;
   [i, pr] = min(ydata(p:p+blockr));
   pr = pr+p-1;
   while (ydata(pr+1) < ydata(pr)) & (pr < length(ydata)), pr = pr+1; end;

   mhi = (2*ydata(p)+ydata(pl)+ydata(pr))/4;
   [i, lm] = min(abs(ydata(pl:p-1)-mhi));
   [i, rm] = min(abs(ydata(p+1:pr)-mhi));
   lm = lm + pl-1;
   rm = rm + p;
   width = xdata(rm)-xdata(lm);
   if ydata(pl)<ydata(pr), bottom = pl; else, bottom = pr; end;
   Q(ct) = xdata(p)/width; 
   mfreq(ct) = xdata(p);
   Pfreq(ct) = sum(ydata(pl:pr))/sum(ydata);   

   if qplot==1;
   hold on
   %plot([xdata(lm) xdata(lm)],[ydata(bottom) ydata(p)],'r: ');
   %plot([xdata(rm) xdata(rm)],[ydata(bottom) ydata(p)],'r: ');
   plot(xdata(p),ydata(p),symbol);
   %set(gca,'fontsize',10);
   %text(xdata(p),ydata(p),'            ');
   %if Pfreq(ct)*100 >= 1,
   %text(xdata(p),ydata(p),[' ',num2str(1/xdata(p),3),' ',num2str(Pfreq(ct)*100,2),'%']);
   %else
   H = text(xdata(p),ydata(p),[' ',num2str(1/xdata(p),3)]);
   set(H,'fontsize',8); 
   %end;
   %plot(xdata(pl),ydata(pl),'b+');
   %plot(xdata(pr),ydata(pr),'b+');
   end;
   clear pry ply i pr pl i lm rm p hi mhi width;
end;  %stop freq for statement.




