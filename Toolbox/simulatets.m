%%simulatets.m
%%generate a time series with a given spectrum by using random number
%%generatino of Fourier coefficients at a fixed resolution
%%assumption is that Phi is given at a coarser resolution from 0 (0 power  )
%% to 1/(2*deltat)
%%Phifreq is provided
%%Use basic notation of Hamming, p. 510. his N=n here.
%%C. Wunsch February 1995. Aims for consistency with convetnions of
%%spectrum2. See note at end.

%%generate a basic frequency scale:
L=(2*n)*deltat;
dfreq=1/L;
freq1=[0:dfreq:n/L];

  %%note must make sure that Phifreq actually gets to the highest
  %%required frequency. E.g. if deltat=1, spectrum2 does not normally
   %%go all the way to 0.5 and need to extrapolate slightly.(if use
   %%a numericallycomputedspectrum. otherwise no problem.
   %%BUT	 the multitaper spectrum does to to 0.5 if use default
   %%option in nfft. otherwise have same problem.
Phiinterp=table1([Phifreq,Phi],freq1);
 %%power density in spectrum2 is sum over M of the an^2
   %% Npts*deltat/M. So mean an^2 is Phi/(npts*deltat);
	
a=randn(n+1,1).*sqrt(Phiinterp/(2*n*deltat)); 
   %%the power is in sine, and half in cosine.
b=randn(n+1,1).*sqrt(Phiinterp/(2*n*deltat));



%%make complex coefficients:
alpha=a-sqrt(-1)*b; %%CHECK THE SIGN CNvnetion
%%make the mean real and zero:
 alpha(1)=0;
%%make the Nyquist frequency real"
  alpha(n)=real(alpha(n));

%%set up to agree with ifft expectations
alpha=n*[alpha;conj(flipud(alpha(2:length(alpha)-1)))];%%factor 
   %of 2 removed because for a unit sine wave fft gives a coefficient
   %%of npts/2 for each of the positive and negative freq. components.
xhat=real(ifft(alpha));

%%%NOTE. set yy=xhat and run spectrum2.

