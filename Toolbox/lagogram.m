%Lagogram.  Lagged cross-correltation at incremental increase times.
%
%     tx: time for x,used for calculating mean time of segements.
%      x: base function.
%   lags: vector of time offsets to compute lagged correlations at
% window: function to window record segements with; also determines
%         segment length.
%
%function [xc,tm]=lagogram(tx,x,ty,y,lags,window);

function [xc,tmx,xvar,yvar]=lagogram(tx,x,ty,y,lags,window,jump);

window=window(:)';
x=x(:);    y=y(:);
tx=tx(:); ty=ty(:);

dt=abs(diff(tx(1:2)));
nwindow=length(window);

%Parse function;
   for ct=1:length(ty)-round(nwindow/dt)-1;
    ypl=ct:ct+round(nwindow/dt)-1;
    yw(ct,:)=(y(ypl)-mean(y(ypl)))';
    yvartemp(ct)=sum(yw(ct,:).^2)/nwindow;
    yw(ct,:)=yw(ct,:).*window;
    yw(ct,:)=yw(ct,:)-mean(yw(ct,:));
    yw(ct,:)=yw(ct,:)/sqrt(sum(yw(ct,:).^2));
   end;

%Compute lagged correlations
ct3=0;
for slide=1:jump:length(ty)-round(nwindow/dt)-1;
    ct3=ct3+1;
    ypl=slide:slide+round(nwindow/dt)-1;
    xpl=ypl+(min(ty)-min(tx));
    tmx(ct3)=mean(tx(xpl));
    xw=x(xpl)-mean(x(xpl));
    xvar(ct3)=sum(xw.^2)/nwindow; 
    xw=xw/sqrt(sum(xw.^2));
    ct2=0;
    for lag=lags;
      ct2=ct2+1;
      if slide-lag>0 & slide-lag<length(yw(:,1)),
      xc(ct3,ct2)=sum(yw(slide-lag,:).*xw');
      yvar(ct3,ct2)=yvartemp(slide-lag);
      else,
        xc(ct3,ct2)=NaN;
      end;    
    end;  
    disp([tmx(ct3)/1000  max(xc(ct3,:))]);
end; %stop the slide loop;

%figure(1); clf; hold on;
%pcolor(tmx,lags,xc');
%shading flat;
%axis tight;
%colorbar vert;
%plot(tmx,zeros(size(tmx)),'k--');
