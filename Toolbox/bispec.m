%Looking at the bispectrum of EOF1. 

if ~exist('V'); load ..\SVD\EOFbp.mat; end;
IAge=13.5:.5:772;

[Bispec, freq] = bicohPH(EOF1,.5,0,'hanning',.1);

pl=find(Bispec<.5);
Bbig=Bispec;
Bbig(pl)=NaN;

figure(1); clf;
contour(freq,freq,Bbig,4);
axis tight;
h=xlabel('Frequency (1/Kyr)'); set(h,'fontsize',13);
h=ylabel('Frequency (1/Kyr)'); set(h,'fontsize',13);
h=title('Bicoherence');         set(h,'fontsize',15);
colorbar('vert');


break;
f  = 1/40.5+[-2:3]*1/98;
fr = round(1000*f)/1000;
set(gca,'ytick',f);
set(gca,'yticklabel',fr);
set(gca,'xtick',f);
set(gca,'xticklabel',fr);
grid on;



