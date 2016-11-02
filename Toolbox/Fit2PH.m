%This functions fits a Cosine wave to a given set of data using fft.
%This one accounts for a red noise background.
%Data must be evenly spaced.  The fit is returned with the same length
%as the data.  Range is the frequency interval to be searched.
%
%function [Fit, amp, freq, phase] = FitPH(Data,dt,range);

function [Fit, amp, freq, phase] = FitPH(Data,dt,range);

time  = 0:dt:(length(Data)-1)*dt;
Data  = Data-mean(Data);
LData = [Data;  zeros(2^15-length(Data),1)];
Lfy   = fft(LData);
Lfy   = Lfy(1:fix(length(LData)/2));
P     = Lfy.*conj(Lfy);
freq  = 1/(length(LData)*dt):1/(length(LData)*dt):length(P)/(length(LData)*dt);
Poly  = polyfit(freq',log(P),5);
Pval  = exp(polyval(Poly,freq));

range_pl = find(freq>=min(range) & freq<=max(range));
freq  = freq(range_pl); P = P(range_pl); 
Lfy   =Lfy(range_pl);   Pval = Pval(range_pl); 
[MP pl]= max(P-Pval');
freq   = freq(pl);
avg    = 7;
phase  = mean(atan2(imag(Lfy(pl-avg:pl+avg)),real(Lfy(pl-avg:pl+avg))));
amp    = 2*sqrt(MP)/length(Data);
Fit    = (amp*cos(freq*time*2*pi+phase))';

%display('Best Fit Cosine: Period, Phase, Amp, std(Data-Fit)');
%display([1/freq phase amp std(Data-Fit)]);

