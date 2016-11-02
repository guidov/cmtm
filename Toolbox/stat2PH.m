%This program takes a record and outputs the Period, Phase, 
%Percent Power, and Q value of the harmonic containing the largest variance.
%This harmonic is then subtracted from the record and the process repeated.
%
%function [Period, Percent, Sig, Q] = stat2PH(x, t, num, frange, qplot);
%
%x    = input data
%t    = time x measuremenst were taken at
%num  = number of iterations to perform.
%frange= frequency space to search over: i.e. 1/800:1/2^14:1/10 or so.
%qplot= do you want these things plotted

function [A, F, phase, x] = stat2PH(x, t, num, frange, qplot);

x=x-mean(x);

%dt =.5;
%t = 2*pi*[0:dt:100];
%t = sort(rand(1,500)*1000);
%x = 2.31*cos(2*pi*t/100+pi*10/180); %+cos(t/41+pi*40/180);
%x = x-mean(x);
%num   = 3;
%qplot = 1; 

for iter=1:num;

way=size(x); if way(1)<way(2); x=x'; end;
way=size(t); if way(1)<way(2); t=t'; end;

if t(1)<0;
[freq,P,H]=fftPH(x,mean(diff(t)),[],14);
pl=find(P==max(P)); pl=pl(1);
A(iter)=2*sqrt(P(pl)*length(x))/length(x);
else;
%freq = 1/800:1/2^11:1/4;
[P, prob, freq]=lombPH(t,x,frange,[]);
pl=find(P==max(P)); pl=pl(1);
A(iter)=sqrt(2*P(pl)/length(x));
end;

F(iter)=freq(pl);

for ct=-179:180;
  xc(ct+180)=sum(cos(2*pi*F(iter)*t+pi*ct/180).*x);
end;  
[y i]=max(xc);
phase(iter)=i-180;
C=A(iter)*cos(2*pi*F(iter)*t+pi*phase(iter)/180);
C=C-mean(C);

if qplot==1;
figure(qplot); clf; subplot(211); hold on;
%fprintf('           Period    Amplitude   Phase\n');
%fprintf('%15.1f %10.1f %10.2f \n',[round(10/F(iter))/10 round(A(iter)*10)/10 phase(iter)]');
plot(freq,P);  
set(gca,'xscale','log'); set(gca,'yscale','log'); axis tight;
subplot(212); hold on;
plot(t,x,'b');
plot(t,C,'r'); axis tight;
title([num2str(1/F(iter),3)]);
drawnow;
end;

x=x-C;
x=x-mean(x);
end;

fprintf('           Period    Amplitude   Phase\n');
fprintf('%15.1f %10.3f %10.0f \n',[1./F' A' phase']');


