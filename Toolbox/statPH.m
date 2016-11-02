%This program takes a power density vector and outputs its statistics
%The statistice are Period, Percent Power, Q value and significance 
%(using ChiPH) for each peak.
%
%function [Period, Percent, Sig, Q] = statPH(freq, P, Ptot, num, v, qplot);
%
%freq = frequencies associated with PSD estimate
%P    = Power Spectral Density (PSD) estimate
%Ptot = Total power in the record.
%num  = number of peaks to look at, sorted by significance.
%v    = degrees of freedom associated with P; for pmtm v = 2*(#windows)-1.
%qplot= do you want these things

function [Period, Percent, Sig, Q] = statPH(freq, P, Ptot, num, v, qplot);


A=size(freq); if A(1)<A(2); freq=freq'; end;
A=size(P);    if A(1)<A(2); P=P'; end;
 

%if qplot>0; figure(qplot); clf; hold on; end;
Ptot=Ptot/std(P);
P=P/std(P);
%Fit an exponential to P
Poly = polyfit(freq,log(P),7);
Pval = exp(polyval(Poly,freq));


%Find max values;
Cmax = [];
for pl =2:length(P)-1; if (P(pl)>=P(pl-1) & P(pl)>P(pl+1));
Cmax = [Cmax pl]; 
end; end;
if Cmax(1)<3; Cmax=Cmax(2:length(Cmax)); end;
if Cmax(length(Cmax))>length(Cmax)-2; Cmax=Cmax(1:length(Cmax)-1); end;

%Find associated mins
Cmin = [];
for ct = 1:length(Cmax); 
  pl = Cmax(ct);  
  while (P(pl-1)<P(pl) & pl>2 & pl<length(P)-2);  pl=pl-1; end;  %To the left;
  Cmin(1,ct)=pl;
  pl = Cmax(ct);
  while (P(pl+1)<P(pl) & pl>2 & pl<length(P)-2);  ;  pl=pl+1; end;  %To the right;
  Cmin(2,ct)=pl;
end;

%Sort according to criteria: Percent Power above noise.
Val = (P(Cmax)-Pval(Cmax))./Pval(Cmax);
[SVal, I] = sort(Val); SVal = flipud(SVal); I = flipud(I);
if length(Cmax)>=num;
SCmax = Cmax(I(1:num));
SCmin = Cmin(:,I(1:num));
else;
  SCmax = Cmax;
  SCmin = Cmin;
end;
Period = 1./freq(SCmax)';
Period = round(Period*10)/10;

%Q value
for ct = 1:length(SCmax);
lmP =mean([P(SCmin(1,ct)) P(SCmax(ct))]);
lf = interp1(P(SCmin(1,ct):SCmax(ct)),freq(SCmin(1,ct):SCmax(ct)),lmP);
rmP =mean([P(SCmax(ct)) P(SCmin(2,ct))]); 
rf = interp1(P(SCmax(ct):SCmin(2,ct)),freq(SCmax(ct):SCmin(2,ct)),rmP);
if qplot>0;
plot(lf,lmP,'g.');
plot(rf,rmP,'g.');
end;
Q(ct) = freq(SCmax(ct))/(rf-lf);
Q     = round(Q*10)/10;
end;

%Percent Power
for ct = 1:length(SCmax);
pl = find(P(SCmin(1,ct):SCmin(2,ct))-Pval(SCmin(1,ct):SCmin(2,ct))>0);
pl = pl+SCmin(1,ct)-1;
Percent(ct) = 100*(sum(P(pl)-Pval(pl)))/Ptot;
Percent = round(Percent*10)/10;
end;

%Significance

for ct = 1:length(SCmax);
Sig(ct) = ChiPH(P(SCmax(ct))/Pval(SCmax(ct)),v);  
end;
Sig = round(100*Sig)/100;

[y i]=sort(Period);
if qplot~=0;
fprintf('           Period    %%Power   Significance   Q\n');
fprintf('&%15.1f& %10.1f& %10.2f& %10.1f\\\\\n',[Period(i)' Percent(i)' Sig(i)' Q(i)']');
end;

if qplot>0;  hold on;
plot(freq,Pval*.041,'r:');
plot(freq,Pval,'k:');
plot(freq,Pval*3.2,'r:');
plot(freq,P);
plot(freq(SCmax),P(SCmax),'r.');
%plot(freq(SCmin),P(SCmin),'k.');
set(gca,'yscale','log'); axis tight;
end;



