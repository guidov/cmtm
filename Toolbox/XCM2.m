Return-Path: <sci14@knorr.whoi.edu>
Received: from pacific-carrier-annex.mit.edu by po9.mit.edu (8.9.2/4.7) id NAA01107; Fri, 25 Jan 2002 13:04:05 -0500 (EST)
Received: from striker.whoi.edu (striker.whoi.edu [128.128.96.129])
	by pacific-carrier-annex.mit.edu (8.9.2/8.9.2) with ESMTP id NAA26929
	for <phuybers@mit.edu>; Fri, 25 Jan 2002 13:03:48 -0500 (EST)
Received: (from knorr@localhost)
	by striker.whoi.edu (8.9.3/8.9.3) id SAA05915
	for phuybers@mit.edu; Fri, 25 Jan 2002 18:03:46 GMT
Received: from localhost (sci14@localhost) by knorr.whoi.edu (8.8.5/8.7.3.scott) with ESMTP id RAA04607 for <phuybers@mit.edu>; Fri, 25 Jan 2002 17:38:19 GMT
Date: Fri, 25 Jan 2002 17:38:19 +0000 (GMT)
From: Peter Huybers <sci14@knorr.whoi.edu>
To: <phuybers@mit.edu>
Subject: xcPHb.m
Message-ID: <Pine.LNX.4.30.0201251737540.4548-100000@mike.knorr.whoi.edu>
MIME-Version: 1.0
Content-Type: TEXT/PLAIN; charset=US-ASCII


%----------------------------------------------------------------
%----------------------FUNCTIONS---------------------------------
%----------------------------------------------------------------

%---This function perfomrs linear interpolation and can handle NaNs---

function [y2]=interpPH(x,y,x2);
for ct=1:length(x2);
pl=find(x<=x2(ct));
pr=find(x>=x2(ct));
if length(pl)==0 | length(pr)==0;
  y2(ct)=NaN;
else,
	  pl=pl(end); pr=pr(1);
  if pl==pr | x(pl)==x(pr),
    y2(ct)=y(pl);
  else;
y2(ct)=(y(pr)-y(pl))/(x(pr)-x(pl)) * (x2(ct)-x(pl)) + y(pl);
end;
end;
end;


%----------Finds the cross-correlation between two vectors-------
function [xc]=xcVal(y1,y2);

pl=find(isnan(y1)==0 & isnan(y2)==0);
y1=y1(pl)-mean(y1(pl));
y2=y2(pl)-mean(y2(pl));
xc=sum(y1.*y2)/sqrt(sum(y1.*y1)*sum(y2.*y2));


%---------Choose Control Points-------------
function [CPpl, CPall]= chooseCP(CPpl,CPall,Mxc);
global xmax xmin dx ymax ymin dy x1 x2 y1 y2 Mx;

dT=round(100*(dx)/(25*length(CPpl)))/100;
T=dT*10;
draw(Mxc,T,0,CPpl);
dist=sqrt( (dy/500)^2+(dx/500)^2);

qgo=0;
count=length(CPpl)-1;
while qgo==0;
count=count+1;
display(['Choose Control Point ',num2str(count)]);
[xt,yt,button]=ginput(1);
display([xt yt]);
if xt<max(Mx)+.1*(max(Mx)-min(Mx)),
[a CPpl(count)]=min( ((Mx-xt)/dx).^2 + ((y2-yt)/dy).^2);
if sum( find( abs(CPpl(1:end-1)-CPpl(end))<dist ) )>0,
pl=find(CPpl(1:end-1)==CPpl(count));
pl=find(CPpl~=CPpl(count));
CPpl=CPpl(pl);
count=count-2;
end;
draw(Mxc,T,0,CPpl);
drawnow;
else, qgo=1;  end;
end;

CPpl(end+1)=length(Mx);
CPpl(end+1)=1;
CPpl=sort(CPpl);
pl=find(diff(CPpl)~=0);
CPpl=CPpl([1 pl+1]);

CPtemp=CPall;
pl=find(isnan(CPall(:,4))==0 | isnan(CPall(:,5)));
CPall=ones(length(Mx(CPpl)),5)*NaN;
CPall(1:length(CPpl),1:3)=[[1:length(CPpl)]', Mx(CPpl)', y2(CPpl)'];

for ct=1:length(pl);
pl2=find(CPtemp(pl(ct),2)==CPall(:,2));
CPall(pl2,4)=CPtemp(pl(ct),4);
CPall(pl2,5)=CPtemp(pl(ct),5);
end;





%-----------Refine Control Points--------------------
function [CPall,cps]=refineCP(CPpl,figNumber,CPall,cps);
global cps xmax xmin dx ymax ymin dy x1 x2 y1 y2 Mx;

heading=['      CP     x-axis     yaxis    left limit   right limit'];
display(heading);
display(CPall);
cpq=1;
while cpq~=0
cpq=input(['Which CP (1-',num2str(length(CPpl)),', 0 to exit, -1 to list
all, -2 set stretch) ']);
if length(cpq)==0, cpq=0; end;
if cpq==0; figure(figNumber); end;
if cpq==-2,
   cps(1)=input('Enter max compaction factor ');
   cps(2)=input('Enter max expansion  factor ');
   cps=abs(cps);
if cps(1)<1 | cps(2)<1; display('you should increase factor(s) which are
less than 1'); end;
end;
if cpq==-1,
display(heading);
display(CPall);
end;
if min(abs([1:length(CPpl)]-cpq))==0,
   display(heading);
   display(CPall(cpq,:));
   cpql=input('Enter left limit  ');
   cpqr=input('Enter right limit ');
   if isempty(cpql), cpql=NaN; end;
   if isempty(cpqr), cpqr=NaN; end;
   CPall(cpq,4)=cpql;
   CPall(cpq,5)=cpqr;
   display(CPall(cpq,:));
end;
end;
for ct=1:length(CPpl);
    if CPall(ct,2)<CPall(ct,4) | CPall(ct,2)>CPall(ct,5),
       CPall(ct,2)=mean(CPall(ct,4:5));
    end;
end;
%stretch core to new length.
Mx=Mx*(max(CPall(:,2))-min(CPall(:,2)))/(max(Mx)-min(Mx));
Mx=Mx-min(Mx)+min(CPall(:,2));



%-------------Maximize cross-correlation---------------
function [Mxc,CPall]=maxXC(CPpl,CPall,figNumber,fixAxis);
global cps xmax xmin dx ymax ymin dy x1 x2 y1 y2 Mx;

ynew=interpPH(Mx,y2,x1);
Mxc=xcVal(ynew,y1);


set(figNumber,'name','Maximizing the Cross Correlation');
display(['Seeking maximum cross-correlation.              ';'Push stop
when you are satisfied with the match.']);

dT=round(100*(dx)/(25*length(CPpl)))/100;
T=dT*10;
if T>min(diff(CPall(:,2))); T=min(diff(CPall(:,2))); end;
[a2Hndl]=draw(Mxc,T,0,CPpl);

count=0;
while get(gca,'Userdata')~=6 & count<5000;
count=count+1;

count2=0;
%while qgo==0;
CPnew=CPall(:,2)+T*randn(size(CPall(:,2)));
for ct=1:length(CPpl);
  if CPnew(ct)<CPall(ct,4) | CPnew(ct)>CPall(ct,5),
     CPnew(ct)=Mx(CPpl(ct));
  end;
end;

qgo=0;
stretch=diff(CPnew')./diff(x2(CPpl));
while sum(find(stretch<1/cps(1) | stretch>cps(2)))>0 & qgo==0;
count2=count2+1;
pl=find(stretch<1/cps(1) | stretch>cps(2));
CPnew([pl pl+1])=CPall([pl pl+1],2)+T*randn(size(CPnew([pl pl+1])));
stretch=diff(CPnew')./diff(x2(CPpl));
if count2>10000;
display(['Cannot find CPs which meet given parameters,     '; ...
         'suggest reducing CPs, changing left/right limits,'; ...
         'or decreasing T.  Probably best to quite program.']);
CPnew=CPall(:,2); qgo==1; count2=0; end;
end;

%end;


xnew =interpPH(x2(CPpl),CPnew',x2);
if sum(isnan(xnew))>0, keyboard; end;

ynew=interpPH(xnew,y2,x1);
xc=xcVal(ynew,y1);
if rem(count,1)==0,
  set(a2Hndl,'String',['attempt ',num2str(count)]);
  drawnow;
  if get(gca,'userdata')==4, T=T+dT;
     draw(Mxc,T,0,CPpl);
     set(gca,'userdata',0);
  end;
  if get(gca,'userdata')==5, T=T-dT;
     if T<dT, T=dT; end;
     draw(Mxc,T,0,CPpl);
     set(gca,'userdata',0);
  end;
end;

if xc>Mxc;
CPall(:,2)=CPnew;
Mx=xnew;
Mxc=xc;
count=0;
[a2Hndl]=draw(Mxc,T,0,CPpl);
drawnow;
end;
end;

set(figNumber,'name','Finished');
display(['Cross-Correlation finished']);

clear *Hndl;

%---------draw plot------------------
function [a2Hndl]=draw(xc,T,attempt,CPpl,drawn);
global cps xmax xmin dx ymax ymin dy x1 x2 y1 y2 Mx;
persistent aHndl rHndl tHndl hp;

xmin=min(Mx);
xmax=max(Mx);
dx=xmax-xmin;


CPpl(end+1)=length(y2);
if nargin<5, delete(hp); end;
hp(1)=plot(x1,y1,'b');
hp(2)=plot(Mx,y2,'r');
hp(3)=plot(Mx,y2,'r.');
hp(4)=plot(Mx(CPpl),y2(CPpl),'ko');
set(hp(4),'markersize',8,'markerfacecolor','k')


if nargin==5;
legend(hp(1:2),'Target','Subject',-1);
aHndl=text('string',['attempt 0'],'position',[xmax+.16*dx
ymax-.12*dy],'fontsize',13,'fontweight','bold');
rHndl=text('string',['r^2
',num2str(round(1000*xc)/1000,3)],'position',[xmax+.16*dx
ymax-.17*dy],'fontsize',10,'fontweight','bold');
tHndl=text('string',['T  ',num2str(T,3)],'position',[xmax+.16*dx
ymax-.22*dy],'fontsize',13,'fontweight','bold');
else,
set(aHndl,'string',['attempt ',num2str(attempt)],'position',[xmax+.16*dx
ymax-.12*dy]);
set(rHndl,'string',['r^2
',num2str(round(1000*xc)/1000,3)],'position',[xmax+.16*dx ymax-.17*dy]);
set(tHndl,'string',['T ',num2str(T,3)],'position',[xmax+.16*dx
ymax-.22*dy]);
end;

a2Hndl=aHndl;


hp(7)=plot([min(Mx) max(Mx)],[ymin-.1*dy ymin-.1*dy],'k');
hp(5)=plot(Mx(CPpl),(ymin-.1*dy)*ones(size(CPpl)),'k.');
stretch=diff(Mx)./diff(x2);
pl1=find(stretch<1);
pl2=find(stretch>=1);
stretch(pl2)=(stretch(pl2)-1)/cps(2);
stretch(pl1) =-(1./stretch(pl1)-1)/cps(1);
hp(6) =plot(Mx(2:end),ymin-.1*dy+.1*dy*(stretch),'k:');
hp(8) =text(min(Mx),ymin-.25*dy,num2str(min(Mx),2));
hp(9) =text(max(Mx)-.06*dx,ymin-.25*dy,num2str(max(Mx),3));
hp(13)=text(mean(Mx),ymin-.25*dy,num2str(mean(Mx),3));
hp(10)=plot([min(Mx) min(Mx)],[ymin ymin-.2*dy],'k');
hp(11)=plot([max(Mx) max(Mx)],[ymin ymin-.2*dy],'k');
hp(12)=plot([mean(Mx) mean(Mx)],[ymin ymin-.2*dy],'k');
hp(14)=text(max(Mx)+.02*dx,ymin,[num2str(cps(2),2),'    Expansion']);
hp(15)=text(max(Mx)+.02*dx,ymin-.2*dy,[num2str(-cps(1),2),'
Compression']);

fixAxis=[xmin-.15*dx xmax+.15*dx ymin-.25*dy ymax+.025*dy];
axis(fixAxis);
drawnow;




