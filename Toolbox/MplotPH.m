%***Modified to use alread interpolated data****************************
%Making a pretty, smoothed Periodogram with Milankovitch lines indicated
%and an estimate of the noise and power at each orbital type.  
%
% function [freq, POx] = MplotPH(Age,Ydata,smooth,name);
%
% POx = smoothed squared Fourier coefficients
% smooth = how many spaces to average POx over, POx is centered.  
% name   = name of core you want to look at, used with getdata.

function [freq, POx] = MplotPH(IAge,IOx,smooth,name);

dt = mean(diff(IAge));
long = length(IAge);
IOx  = IOx'; %';
%IAge = Age(1):dt:Age(length(Age));
%IOx = interp1(Age,Ox,IAge);

%Make the Periodorgram
[freq, POx, H] = fftPH(IOx,dt,[]);
if smooth>1; SPOx = smoothPH(POx,smooth); else; SPOx = POx; end;
H=plot(freq,SPOx,'k'); 
set(H,'linewidth',1.7);
set(gca,'yscale','log');
lplot=find(freq>=1/500 & freq <= 1/22);
axis([1/700 1/22 min(SPOx(lplot)) 1.75*max(SPOx(lplot))]);

%Finding Power at the Milankovitch Freqs and some plotting
Ilines = [1/41.09, 1/39.72, 1/40.39, 1/53.86, 1/28.91;   %Obliquity 
          1/23.71, 1/22.39, 1/22.39, 1/22.39, 1/22.39;   %1/18.96, 1/19.12, 1/16.91; Precession
          1/94.78, 1/123.8, 1/98.72, 1/130.6, 1/404.2];  %Eccentricity 1/404.2 not included

for row=1:3
Irow = Ilines(row,:);
[Nfreq, Pband] = qPHc(POx,freq,3,Irow);
if row==1, OPow = sum(Pband); Opl = find(isnan(Nfreq)==1); end;
if row==2, PPow = sum(Pband); Ppl = find(isnan(Nfreq)==1); end;
if row==3, EPow = sum(Pband); Epl = find(isnan(Nfreq)==1); end;

for col=1:5
Iline = Ilines(row,col);
label = [];
if row==1, color='r:'; 
if col==1 | col==4 | col==5, label = 'Ob';
end; end;  
if row==2, color='g:'; 
if col==1 | col==2 | col==3 | col==5, label = 'Pr';
end; end;
if row==3, color='b:'; 
if col==1 | col==2, label='Ec'; 
end; end;
plot([Iline Iline],[min(POx) 1.05*max(POx)],color); 
if length(label)>0, H=text(Iline,1.35*max(POx),[label,num2str(1/Iline,3)]); set(H,'fontsize',7); end;
end; end;

%Fit a line to the Power coefficients using non-Milankovitch freqs
Nfreq([Opl' Ppl' Epl']) =  NaN; %'
Pplace   = find(isnan(Nfreq)==0);
[Poly S] = polyfit(freq(Pplace),log(POx(Pplace)),3);
Pval     = polyval(Poly,freq);
Pstd     = std(POx(Pplace));
plot(freq,exp(Pval),'k:');
plot(freq,exp(Pval)+2*Pstd,'k:');

%Finding Total Power
PPow = PPow-sum(Pval(Ppl));
EPow = EPow-sum(Pval(Epl));
OPow = OPow-sum(Pval(Opl));

SPow = PPow+EPow+OPow;
NPow = sum(POx)-SPow;
display(['Signal to Noise Ratio: ',num2str(SPow/NPow,3)]);

%Labelling the plot with this info.
h=title(['Power Spectrum of d18O in Core ',name]);
set(h,'fontsize',15)
h=ylabel('Power Density');
set(h,'fontsize',12)
h=xlabel(['Frequency in Cycles/Kyr:  var(d18O)=',num2str(var(IOx),2),' Noise=',num2str(100*NPow/sum(POx),2),'% Obliquity=',num2str(100*OPow/sum(POx),2),'% Precession=',num2str(100*PPow/sum(POx),2),'% Eccentricity=',num2str(100*EPow/sum(POx),2),'%']);
set(h,'fontsize',12)




