%Looking at variance explained by a band.


dt=.1;
x=randn(1000,1);

var=sum(x.^2);

fx=fft(x);   

abs(fx).^2


[freq,P]=fftPH(x,dt,[],0);

pl=find(freq>1/10 & freq<1);  

figure(1); clf; hold on;
plot(freq,P); logPH;
plot(freq(pl),P(pl),'r');

E=sum(P(pl))/sum(P),

