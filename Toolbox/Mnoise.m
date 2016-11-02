%This function make the appropriate noise background:
%This can be extended to any type of noise.
%
%function [Noise] = Mnoise(time,color,qplot);
%
%time = Age and resolution of core
%color = type of noise: 'r'=red, 'w'=white.
%qplot= 0 for no plotting, >0 selectes figure #.


function [Noise] = Mnoise(time,color,qplot);

%color = 'r';
%color = 'w';

%----------Reset Random # Generator---------------%
rand('state',sum(100*clock));
randn('state',sum(100*clock));



%----------Make Frequency Attributes--------------%
dt   = mean(diff(time));
long = length(time);
df   = 1/(long*dt);
freq = [df:df:(df*fix(long/2))];

if color=='w';
A   = 1+randn(1,length(freq));
end;

if color=='r';
pause;
Phfit =  [-6.1181   -8.1759];
Plfit =  [545.1180  -12.9007]; %Constants determined from Lred.m are Plfit and Phfit - two lines.
pll  = find(freq<1/110);
plh  = find(freq>1/110);
Aval = exp([polyval(Plfit,freq(pll)) polyval(Phfit,freq(plh))]);
A    = Aval./freq;
A    = A/std(A);
A    = A.*(1+.25*randn(1,length(A)));
end;


%----------Convert to Time (using random phasing)---%
Noise = sin(freq(1)*2*pi*time+rand(1,1)*2*pi);
for ct = 2:length(freq);
Noise = Noise+A(ct)*sin(freq(ct)*2*pi*time+rand(1,1)*2*pi);
end;

Noise = Noise-mean(Noise);
Noise = Noise/std(Noise);

%----------Plotting (Optional)----------------------%
if qplot>0; 
figure(qplot); clf; subplot(311); 
semilogy(freq,A.^2);
axis tight;
title('Frequency vs. Amplitude^2');
ylabel('Amp^2 (log)');
subplot(313);
plot(time,Noise);
title('Background Noise');
ylabel('Noise (linear)');
xlabel('Time (Kyrs)');
axis tight;
subplot(312);
[F,P,H] = fftPH(Noise,dt,'b');
%pl = find(F>1/500 & F<1/15);
%axis([min(F(pl)) max(F(pl)) min(P(pl)) max(P(pl))]);
axis tight;
title('Periodogram of Noise Signal (Low Freqs Only)');
ylabel('Normalized Amp^2 (log)');
end;
