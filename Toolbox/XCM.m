%Maximizes the cross-correlation between two vectors.
%
%x1 = fixed x ordinate (blue)
%y1 = fixed y oridnate (blue)
%x2 = mallable x ordinate (red)
%y2 = values associated with x2 (red)
%NCP= number of control points
%
%function [MCP, Mx]=XM(x1,y1,x2,y2);

function [MCP, Mx]=XCM(x1,y1,x2,y2);
warning off;
global cps xmax xmin dx ymax ymin dy x1 x2 y1 y2 Mx;
warning on;

[a b]=size(x1);
if a>b, x1=x1'; end;
[a b]=size(y1);
if a>b, y1=y1'; end;
[a b]=size(x2);
if a>b, x2=x2'; end;
[a b]=size(y2);
if a>b, y2=y2'; end;

if nargin~=4 & exist('XCtemp.mat')==2,
  load XCtemp.mat;
else,
save XCtemp.mat;
end;

format compact;
%x2=sort(21*rand(1,101));
%y2=sin([0:.2:20])+rand(size(x2));
%x1=0:.2:20;
%y1=sin([0:.2:20])+rand(size(x2));


xmax=max(x2);
xmin=min(x2);
dx=xmax-xmin;
ymax=max([max(y2) max(y1)]);
ymin=min([min(y2) min(y1)]);
dy=ymax-ymin;

fixAxis=[xmin-.15*dx xmax+.15*dx ymin-.25*dy ymax+.025*dy];

   figNumber=figure( ...
      'Name',['Choose Control Points'], ...
      'NumberTitle','off', ...
      'DoubleBuffer','on', ...
      'Visible','on', ...
      'Color','white', ...
      'BackingStore','off');
   axes( ...
      'Position',[0.05 0.05 0.75 0.90], ...
       'Visible','on',    ...
       'DrawMode','fast', ...
       'Units','normalized', ...
       'NextPlot','add');

   %===================================
   % Information for all buttons
   labelColor=[0.8 0.8 0.8];
   yInitPos=0.90;
   xPos=0.85;
   btnLen=0.1;
   btnWid=0.085;
   % Spacing between the button and the next command's label
   spacing=0.04;

   %====================================   % The CONSOLE frame
   frmBorder=0.02;
   yPos=0.05-frmBorder;
   frmPos=[xPos-frmBorder yPos btnLen+2*frmBorder 0.9+2*frmBorder];
   h=uicontrol( ...
      'Style','frame', ...
      'Units','normalized', ...
      'Position',frmPos, ...
      'BackgroundColor',[0.50 0.50 0.50]);

   %====================================
   % The STOP button
   btnNumber=6;
   yPos=0.90-(btnNumber-1)*(btnWid+spacing);
   labelStr='Stop';
   % Setting userdata to -1 (=stop) will stop the demo.
   callbackStr='set(gca,''Userdata'',6)';

   % Generic button information
   btnPos=[xPos yPos-spacing btnLen btnWid];
   stopHndl=uicontrol( ...
      'Style','pushbutton', ...
      'Units','normalized', ...
      'Position',btnPos, ...
      'Enable','off', ...
      'String',labelStr, ...
      'Callback',callbackStr);


   %====================================
   % The xc button
   btnNumber=1;
   yPos=0.90-(btnNumber-1)*(btnWid+spacing);
   labelStr='XC';
   callbackStr='set(gca,''Userdata'',1)';

   % Generic button information
   btnPos=[xPos yPos-spacing btnLen btnWid];
   xcHndl=uicontrol( ...
      'Style','pushbutton', ...
      'Units','normalized', ...
      'Position',btnPos, ...
      'Enable','off', ...
      'String',labelStr, ...
      'Callback',callbackStr);


   %====================================
   % The CP button (choose control points)
   btnNumber=2;
   yPos=0.90-(btnNumber-1)*(btnWid+spacing);
   labelStr='CP';
   callbackStr='set(gca,''Userdata'',2)';

   % Generic button information
   btnPos=[xPos yPos-spacing btnLen btnWid];
   cpHndl=uicontrol( ...
      'Style','pushbutton', ...
      'Units','normalized', ...
      'Position',btnPos, ...
      'Enable','off', ...
      'String',labelStr, ...
      'Callback',callbackStr);


   %====================================
   % The Refine button (refine control point specs)
   btnNumber=3;
   yPos=0.90-(btnNumber-1)*(btnWid+spacing);
   labelStr='Refine';
   callbackStr='set(gca,''Userdata'',3)';

   % Generic button information
   btnPos=[xPos yPos-spacing btnLen btnWid];
   refHndl=uicontrol( ...
      'Style','pushbutton', ...
      'Units','normalized', ...
      'Position',btnPos, ...
      'Enable','off', ...
      'String',labelStr, ...
      'Callback',callbackStr);

   %====================================
   % The Increase Temperature button
   btnNumber=4;
   yPos=0.90-(btnNumber-1)*(btnWid+spacing);
   labelStr='+T';
   callbackStr='set(gca,''Userdata'',4)';

   % Generic button information
   btnPos=[xPos yPos-spacing btnLen btnWid];
   TplusHndl=uicontrol( ...
      'Style','pushbutton', ...
      'Units','normalized', ...
      'Position',btnPos, ...
      'Enable','off', ...
      'String',labelStr, ...
      'Callback',callbackStr);


   %====================================
   % The Decrease Temperature button
   btnNumber=5;
   yPos=0.90-(btnNumber-1)*(btnWid+spacing);
   labelStr='-T';
   callbackStr='set(gca,''Userdata'',5)';

   % Generic button information
   btnPos=[xPos yPos-spacing btnLen btnWid];
   TminHndl=uicontrol( ...
      'Style','pushbutton', ...
      'Units','normalized', ...
      'Position',btnPos, ...
      'Enable','off', ...
      'String',labelStr, ...
      'Callback',callbackStr);



   %====================================
   % The CLOSE button
   labelStr='Close';
   callbackStr='set(gca,''userdata'',7);';
   closeHndl=uicontrol( ...
      'Style','push', ...
      'Units','normalized', ...
      'Position',[xPos 0.05 btnLen 0.10], ...
      'String',labelStr, ...
      'Callback',callbackStr);

   % Uncover the figure
   hndlList=[stopHndl closeHndl];
   set(figNumber,'Visible','on', ...
      'UserData',hndlList);

   cla;
   axHndl=gca;
   figNumber=gcf;
   hndlList=get(figNumber,'Userdata');
   stopHndl=hndlList(1);
   closeHndl=hndlList(2);
   set(closeHndl,'Enable','off');
   set(stopHndl,'Enable','off');
   set(axHndl, ...
      'UserData',0, ...
      'DrawMode','fast', ...
      'Visible','off');


% ================ Start the Cross Correlation Program==================== %

ynew=interpPH(x2,y2,x1);
pln=find(isnan(ynew)==0);
xc=xcVal(ynew(pln),y1(pln));
cps=[2 2];
Mxc=xc;
%keyboard;
Mx=x2;
qb=0;
q=0;
qgo=0;
CPpl=[1 length(x2)];
CPall=[[1; 2] x2(CPpl)' y2(CPpl)' [NaN; NaN] [NaN; NaN]];
draw(Mxc,0,0,CPpl,0);
display(['Push XC to begin correlation             '; ...
         'Push CP to select control points         '; ...
         'Push Ref to refine control point behavior']);
axis(fixAxis);
set(gca,'userdata',0);
set([closeHndl cpHndl xcHndl refHndl],'Enable','on');
set([TplusHndl TminHndl stopHndl],'Enable','off');

while q~=7;
  if q==1; set([stopHndl TplusHndl TminHndl],'Enable','on');
           set([refHndl xcHndl cpHndl closeHndl],'Enable','off');
           [Mxc,CPall]=maxXC(CPpl,CPall,figNumber,fixAxis);
           set([refHndl cpHndl xcHndl closeHndl],'Enable','on');
           set([stopHndl TplusHndl TminHndl],'Enable','off');
           set(gca,'userdata',0);
           end;
  if q==2; [CPpl CPall]=chooseCP(CPpl,CPall,Mxc);
           set(gca,'userdata',0);  q=0;
           end;
  if q==3; CPall=refineCP(CPpl,figNumber,CPall,cps);
           set(gca,'userdata',0); q=0;
           end;
drawnow;
q=get(gca,'userdata');
end;

display('Saving file to XCout.txt and XCout.mat');
Data=ones(max([length(x1) length(x2)]),5)*NaN;
Data(:,1:2)=[x1' y1'];
Data(1:length(x2),3:5)=[x2' Mx' y2']; 
ControlPoints=[NaN*ones(length(CPall(:,1)),1) CPall];
save XCout.txt ControlPoints Data -ASCII  
keyboard;
CPLabel=['Control Point #  x2     y2          left limit    right limit'];
DataLabel=['x1  y1 x2  x2new  y2'];
save XCout.mat Data CPall CPLabel DataLabel;
clear global;
display('Finished');
close(gcf);


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
cpq=input(['Which CP (1-',num2str(length(CPpl)),', 0 to exit, -1 to list all, -2 set stretch) ']);
if length(cpq)==0, cpq=0; end;
if cpq==0; figure(figNumber); end;
if cpq==-2,
   cps(1)=input('Enter max compaction factor ');
   cps(2)=input('Enter max expansion  factor ');
   cps=abs(cps);
if cps(1)<1 | cps(2)<1; display('you should increase factor(s) which are less than 1'); end;
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
display(['Seeking maximum cross-correlation.              ';'Push stop when you are satisfied with the match.']);

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
rHndl=text('string',['r ',num2str(round(1000*xc)/1000,3)],'position',[xmax+.16*dx ymax-.17*dy],'fontsize',10,'fontweight','bold');
tHndl=text('string',['T  ',num2str(T,3)],'position',[xmax+.16*dx
ymax-.22*dy],'fontsize',13,'fontweight','bold');
else,
set(aHndl,'string',['attempt ',num2str(attempt)],'position',[xmax+.16*dx
ymax-.12*dy]);
set(rHndl,'string',['r ',num2str(round(1000*xc)/1000,3)],'position',[xmax+.16*dx ymax-.17*dy]);
set(tHndl,'string',['T ',num2str(T,3)],'position',[xmax+.16*dx ymax-.22*dy]);
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
hp(15)=text(max(Mx)+.02*dx,ymin-.2*dy,[num2str(-cps(1),2),'Compression']);

fixAxis=[xmin-.15*dx xmax+.15*dx ymin-.25*dy ymax+.025*dy];
axis(fixAxis);
drawnow;



