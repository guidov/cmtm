%function [s, c, ph, ci, phi] = cmtm(x,y,dt,NW,qplot,confn,qdet);
%
%Multi-taper method coherence estimate using an even weighting for the 
%tapered estimates.  In addition, a Monte Carlo estimation procedure is
%available to estimate coherence and phase 95% confidence limits.  An effort
%is made to correct for the bias ineherent to coherence estimates (see 
%discussion on line 126).
%
% Inputs:
%      x       - Input data vector 1.
%      y       - Input data vector 2.
%      dt      - Sampling interval (default 1)
%      NW      - Number of windows to use (default 8)
%      qplot   - Plot the results, yes=1, No=0 (default).  The upper tickmarks
%                  indicate the bandwidth of the coherence and phase estimates.                  
%      confn   - Number of Monte Carlo iterations to use (default 50)
%                  in computing approximate 95% confidence intervals. 0 or
%                  omit to forgo confidence intervals.  In some cases as many
%                  as 200 iterations are required for results to converge.
%                  Stochastic realizations assume the spectral power of x and y 
%                  can be approximated as a second order polynomial.  
%      qdet    - To retain the original time-series x (e.g. if x is
%                  deterministic) and only form a stochastic realization  
%                  of y use qdet=1, otherwise qdet=0 (default).
%
% Outputs:
%      s       - frequency
%      c       - coherence
%      ph      - phase
%      ci      - 95% coherence confidence level
%      phi     - 95% phase confidence interval, bias corrected
%
%
%Peter Huybers
%MIT, 2003
%phuyber@mit.edu

function [s, c, ph, ci, phi] = cmtm(x,y,dt,NW,qplot,confn,qdet,fnum);

x   = x(:)-mean(x); y=y(:)-mean(y);

%check input
if nargin<7, qdet=0;   end;
if nargin<6, qplot=0;  end;
if nargin<5, confn=0;  end; 
if length(confn)==0, confn=0; end;
if nargin<4, NW=8;     end;
if length(NW)==0, NW=8;end;
if nargin<3, dt=1;     end;
if length(dt)==0, dt=1;end;

%define some parameters
N   = length(x);
k   = min(round(2*NW),N); 
k   = max(k-1,1);
s   = (0:1/(N*dt):1/dt-1/(N*dt))';
pls=2:(N+1)/2+1;
if rem(length(y),2)==1; pls=pls(1:end-1); end;

%Compute the discrete prolate spheroidal sequences, requires the spectral analysis toolbox.
[E,V]=dpss(N,NW,k);

%Compute the windowed DFTs.
fkx=fft(E(:,1:k).*x(:,ones(1,k)),N);
fky=fft(E(:,1:k).*y(:,ones(1,k)),N);

if confn>0, 
  Fx=mean(abs(fkx)')';
  Fy=mean(abs(fky)')';
end;

%Normalize the energy of each tapered transform
for ct=1:length(E(1,:));
  fkx(:,ct)=fkx(:,ct)/std(fkx(:,ct));
  fky(:,ct)=fky(:,ct)/std(fky(:,ct));
end;

%Compute coherence
Cxy= sum( (fkx.*conj(fky))' );
c  = abs(Cxy)./sqrt(sum(abs(fkx').^2).*sum(abs(fky').^2));
ph = angle(Cxy)*180/pi;

%Estimate confidence intervals
if confn>0,
  disp('Estimating confidence intervals');
  for iter=1:confn;
    if rem(iter,10)==0, disp(['iteration: ',num2str(iter)]); end;
    
    %Coherence uncertainty estimates
    if qdet==0,
      fxc=abs(fft(x));
      pl=find(fxc(pls)<2.5*std(fxc))+1;
      fxc2=exp(polyval(polyfit(log(s(pl)),log(fxc(pl)),2),log(s(pls))));
      if rem(length(x),2)==1;
	fxc2=[0; fxc2; flipud(fxc2)];
      else,
	fxc2=[0; fxc2; flipud(fxc2(2:end))];
      end;	
    end;  
    
    fyc=abs(fft(y));
    pl=find(fyc(pls)<2.5*std(fyc))+1;
    fyc2=exp(polyval(polyfit(log(s(pl)),log(fyc(pl)),2),log(s(pls))));
    if rem(length(y),2)==1;
      fyc2=[0; fyc2; flipud(fyc2)];
    else,
      fyc2=[0; fyc2; flipud(fyc2(2:end))];
    end;	
    ys=randn(size(y)); 
    ys=ys-mean(ys); ys=ys/std(ys);
    ys=real(ifft(fft(ys).*fyc2));    
    if qdet==0
      xs=randn(size(x));
      xs =real(ifft(fft(xs).*fxc2));
    else,
      xs=x;
    end;
    [si, ci(iter,:), dum]=cmtm(xs,ys,dt,NW);

    %Phase uncertainty
    fx=fft(xs);
    fx=fx/sum(abs(fx));
    fy=fft(randn(size(y))+rand).*fft(y); 
    fy=fy/sum(abs(fy));
    cb=c-(1-c).^3; 
    %NOTE: coherence is a biased estimator. Runs with known signals indicated
    %a positive bias of roughly +.3 for incoherent processes, with the biass
    %tapering off toward higher true coherence.  However, results vary
    %according to the number of windows used and the spectral structure of the
    %record.  An empirically derived adjustment is made to correct for the
    %bias.  This adjustment appears to be approximately correct over a range of
    %signals and numbers of windows.    
    pl=find(cb<0); 
    cb(pl)=0; 
    ys =real(ifft(fy.*sqrt(1-cb'.^2)));    
    ys =ys+real(ifft(fx.*cb')); 
    [si, ciph(iter,:), phi(iter,:)]=cmtm(xs,ys,dt,NW);
  end;

  %sorting and averaging to determine confidence levels.
  pl=round(.95*iter);  
  ci=sort(ci);    
  ci=ci(pl,:);
  if qdet==0,
    ci=mean(ci)*ones(size(ci));
  end;
  pl=round(.975*iter);  
  phi=sort(phi);  
  phi=[phi(pl,:); -phi(iter-pl+1,:)];
  phi=mean(phi);
  %phase uncertainty is smoothed with a triangle window, 
  %except at the edges where boxcar averaging is used.
  phi(1:3)=mean(phi(1:3));
  phi(end-2:end)=mean(phi(end-2:end));
  temp=conv(phi,[.25 .5 .25]);
  phi(4:end-3)=temp(5:end-4);
  phi=[phi; -phi]; 
  phi=phi+repmat(ph(pls),2,1);
end;

%Cut to one-sided funtions
c = c(pls);
s = s(pls)';
ph=ph(pls);

%plotting
if qplot==1,
  %coherence
  %figure(gcf); clf;
  subplot(2,2,fnum); hold on;
  plot(s,c);
  h=ylabel('coherence'); 
  if confn>0;
    plot(si,ci,'k--');
    pl=find(c>ci);
    title([num2str(100*length(pl)/length(c),2),'% of estimates above 95% confidence level']);
  end;
  axis tight; h=axis; axis([h(1:2) 0 1.025]); h=axis;
  w  = NW/(dt*N);   %half-bandwidth of the dpss
  plot([s(1) h(2)],[1.02 1.02],'k');
  for ds=min(s):2*w:max(s);
    plot([ds ds],[.98 1.02],'k');
  end;
  h=axis; axis([h(1:3) 1.025]);
  set(gca,'xscale','log');
  
  %phase
  subplot(2,2,fnum+1); hold on;
  plot(s,ph);
  if confn>0,
    cu=phi(1,:);
    cl=phi(2,:);
    col=[.9 .9 .9];
    h=fill([s(1) s(1:end) fliplr([s(1:end) s(end)])],[cu(1) cl fliplr([cu cl(end)])],col);
    set(h,'edgecolor',col);

    pl=find(cu<=180); cu(pl)=-180;
    pl=find(cu> 180); cu(pl)=cu(pl)-360;
    clt=-180*ones(size(cl));
    h=fill([s(1) s(1:end) fliplr([s(1:end) s(end)])],[cu(1) clt fliplr([cu clt(end)])],col);
    set(h,'edgecolor',col);

    pl=find(cl>=-180); cl(pl)=180;
    pl=find(cl< -180); cl(pl)=cl(pl)+360;
    cut=180*ones(size(cl));
    h=fill([s(1) s(1:end) fliplr([s(1:end) s(end)])],[cut(1) cl fliplr([cut cl(end)])],col);    
    set(h,'edgecolor',col);        
  end;
  h=plot(s,ph); %set(h,'linewidth',1.5);      
  plot(s,zeros(size(s)),'k--');
  axis tight; h=axis; axis([h(1:2) -180 180]);
  set(gca,'xscale','log');
  xlabel('frequency (cycles/deltat)')
  h=ylabel('phase');
end;









