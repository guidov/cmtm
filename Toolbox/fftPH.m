%A periodogram using the fft function
%
%function [freq,P,H] = fftPH(ydata,deltat,linetype,pad)
%
%Pad      = choose what power of two you would like to pad to (default is 13).
%H        = handle to the graphics that are displayed
%linetype = color and linestyle for plot e.g. 'r:' for red, dotted.
%           Leave empty, [], for no plot.

function [freq,P,H,Ang] = fftPH(ydata,deltat,linetype,pad)

if length(ydata(1,:))==1; ydata=ydata'; end;
long = length(ydata);
if nargin<3, linetype=[]; end;
if nargin<4, pad=0; end;

ydata = ydata-mean(ydata);
Lydata = [ydata  zeros(1,2^pad-length(ydata))];
fy = fft(Lydata);
fy = fy(2:length(fy));
fy = fy(1:fix(length(fy)/2));
P  = fy.*conj(fy); 
Ang= angle(fy)*180/pi; 
df = 1/(length(Lydata)*deltat);
freq = df:df:df*length(P);
P = P*sum(ydata.^2)/(sum(P));
P=P/length(ydata);
%keyboard;

%Only use lowest meaningfull freq.
%pl=find(freq>=1/(length(ydata)*deltat));
%P=P(pl); freq=freq(pl);
if ~isempty(linetype)
   if length(linetype)>3; 
     pl=str2num(linetype(4:6)); ph=str2num(linetype(7:9));
     pl=find(freq<1/ph & freq>1/pl); P=P(pl); freq=freq(pl); 
   end;
   hold on;
   H = plot(freq,P,linetype(1));
   h=title(['Periodogram ', date]);  set(h,'fontsize',14); 
   h=ylabel('Power (Units^2/Cycle/dt)');  set(h,'fontsize',12); 
   h=xlabel('Frequency'); set(h,'fontsize',12); 
   axis tight; h=axis;

   %Make the 95% confidence interval
   ci=chi2confPH(.95,2);
   Poly = polyfit(freq,log(P),3);
   Pval = exp(polyval(Poly,freq));  
   plot(freq,Pval*ci(1),'r:'); 
   plot(freq,Pval*ci(2),'r:'); 
   plot(freq,Pval,'r'); 
else
   H = [];
end;
freq = freq;   
