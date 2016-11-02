%This function looks for a central peak and then calculates the period of an oscillation 
%in sedimentation that would be required to account for a shift from the central peak to 
%the observed neighbors
%
%[FP Per] = modPH(freq,P,CP);
%
%FP = the frequency that the surrounding peaks are located at
%Per= the periodicity required to account for the presumed shift
%P  = the power from a power spectrum
%CP = the central peak to adjust from 

clear; clc; format compact;
figure(1); clf; hold on;


CP = 40.5;
F = 1/41.1;
F2= 1/23.8;       %for ODP980 900 to 2200 with CP = 40.5 and AOx and CF as defined
MaxTime = 2100;   %for ODP677 900 to 2200 works            
MinTime = 850;    %for ODP849 850 to 2100 works     
name    = 'ODP849';
[Dep Age Ox Car] = getdata(name);

p1    = find(Age<=MaxTime & Age>=MinTime);
Age   = Age(p1); Dep = Dep(p1); Ox = -Ox(p1); clear Car;
AOx   = detunePH(Age,Dep,Age(length(Age)));  clear Age;
%AOx   = AOx*40.3/41.7; %used for ODP677
%AOx    = AOx*40.3/41.2;  %used for ODP849
%AOx   = AOx*40.3/38.6;  %used for ODP980
%AOx   = AOx*41.1/40.1;  %used for ODP677
dt    = mean(diff(AOx));
IAge  = min(AOx):mean(diff(AOx)):max(AOx);
IOx   = interp1(AOx,Ox,IAge);
ID    = interp1(AOx,Dep,IAge);
[freq P] = fftPH(diff([0 IOx]),dt,[]);
place = find(freq > 1/800 & freq<1/13);
P = P(place);
freq = freq(place);
P = smoothPH(P,5);
subplot(1.5,1,1); hold on;
H=plot(freq,P,'b')
set(gca,'yscale','log');
axis tight;

%find top/bottom in filter 18O data
Pmax = [];
for pl =2:length(P)-1;
if P(pl)>=P(pl-1) & P(pl)>P(pl+1), Pmax=[Pmax; pl]; end;
end; 

[a Cpl] = min(abs(1/CP-freq(Pmax)));
CF = freq(Pmax(Cpl));
CF =  1/40.45;  %for ODP980 and DSDP607
plot(freq(Pmax),P(Pmax),'b.');
plot(CF,P(Pmax(Cpl)),'r*');

DP = 1./(abs(freq(Pmax)-CF));
for ct=1:length(DP); if ct~=Cpl;
H=text(freq(Pmax(ct)),P(Pmax(ct)),[' ',num2str(DP(ct),3)]);
set(H,'fontsize',8);
end; end;
H=text(CF,P(Pmax(Cpl)),[' ',num2str(1/CF,3)]);
set(H,'fontsize',12);

%Finding Power at the Milankovitch Freqs and some plotting
Ilines = [1/41.09, 1/39.72, 1/40.39, 1/53.86, 1/28.91;   %Obliquity 
          1/23.71, 1/22.39, 1/22.39, 1/22.39, 1/22.39;   %1/18.96, 1/19.12, 1/16.91; Precession
          1/94.78, 1/123.8, 1/98.72, 1/130.6, 1/404.2];  %Eccentricity 1/404.2 not included

for row=1:3; Irow = Ilines(row,:);
for col=1:5; Iline = Ilines(row,col);
label = [];
if row==1, color='r:'; if col==1 | col==4 | col==5, label = 'Ob';
end; end;  
if row==2, color='g:'; if col==1 | col==2 | col==3 | col==5, label = 'Pr';
end; end;
if row==3, color='b:'; if col==1 | col==2, label='Ec'; 
end; end;
plot([Iline Iline],[min(P) 1.05*max(P)],color); 
if length(label)>0, H=text(Iline,1.35*max(P),[label,num2str(1/Iline,3)]); set(H,'fontsize',7); end; end; end;




subplot(4,1,4);
plot(IAge,IOx);
xlabel('Kyrs');
title(['Modulation of Core ',name]);
axis tight
