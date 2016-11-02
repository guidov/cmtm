%function [s, c, ph, ci, phi] = cmtmPH(x,y,dt,NW,confn,qplot);
%
%simple multi-taper method coherence estimate and Monte Carlo
%confidence estimation procedure.
%
% Inputs:
%      x        - Input data vector 1.
%      y        - Input data vector 2.
%      dt       - sampling interval
%      NW       - number of windows to use (8 is often a good choice)
%      confn    - number of Monte Carlo iterations to use (default 50)
%                 in computing approximate 95\% confidence intervals. 0 or
%                 omit to forgo confidence intervals.  In some cases as many
%                 as 200 iterations are required for results to converge.
%      qplot    - plot the results, 1= Yes, 0 or omit = No.  The upper
%                 tickmarks indicate the bandwidth of the coherence and
%                 phase estimates.
%
% Outputs:
%      s       - frequency
%      c       - coherence
%      ph      - phase
%      ci      - 95% coherence confidence level
%      phi     - 95% phase confidence interval
%
%
%Peter Huybers
%MIT, 2003
%phuyber@mit.edu

function [s, c, ph, ci, phi] = cmtmPH(x,y,dt,NW,confn,qplot);

x   = x(:)-mean(x); y=y(:)-mean(y);

%check input
if nargin<6, qplot=0;  end;
if nargin<5, confn=0; end; 
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

%Compute the discrete prolate spheroidal sequences, requires the
%spectral analysis toolbox.
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
c  = abs(Cxy)./sqrt(sum([fkx.*conj(fkx)]').*sum([fky.*conj(fky)]'));
ph = angle(Cxy)*180/pi;

%Estimate confidence intervals
if confn>0,
  disp('Estimating confidence intervals');
  for iter=1:confn;
    if rem(iter,10)==0, disp(['iteration: ',num2str(iter)]); end;

    %Make filters for phase uncertainty estimates
    fxc=abs(fft(x));
    pl=find(fxc(pls)<2.5*std(fxc))+1;
    fxc2=exp(polyval(polyfit(log(s(pl)),log(fxc(pl)),2),log(s(pls))));
    if rem(length(x),2)==1;
      fxc2=[0; fxc2; flipud(fxc2)];
    else,
      fxc2=[0; fxc2; flipud(fxc2(2:end))];
    end;	
      
    fyc=abs(fft(y));
    pl=find(fyc(pls)<2.5*std(fyc))+1;
    fyc2=exp(polyval(polyfit(log(s(pl)),log(fyc(pl)),2),log(s(pls))));
    if rem(length(y),2)==1;
      fyc2=[0; fyc2; flipud(fyc2)];
    else,
      fyc2=[0; fyc2; flipud(fyc2(2:end))];
    end;	

    %coherence
    ys=randn(size(y)); 
    ys=ys-mean(ys); ys=ys/std(ys);
    ys=real(ifft(fft(ys).*fyc2));    
    xs=randn(size(x)); 
    xs=xs-mean(xs); xs=xs/std(xs);
    xs =real(ifft(fft(xs).*fxc2));
    [si, ci(iter,:), dum]=cmtmPH(xs,ys,dt,NW);

    %phase
    fy=fft(randn(size(y))+rand).*fft(y); 
    fx=fft(randn(size(x))+rand).*fft(x); 
    fy=fy/sum(abs(fy));
    fx=fx/sum(abs(fx));
    cb=c-(1-c).^3; pl=find(cb<0); cb(pl)=0; 
    ys =real( ifft(fy.*sqrt(1-cb'.^2)));    
    ys =ys+real(ifft(fx.*cb'));    %NOTE: coherence is a biassed estimator.
                                   %Runs with known signals indicated a bias
                                   %of +.3 for incoherent processes, with
                                   %the biass tapering off toward higher
                                   %true coherence.  Thus an empiracally
                                   %derived  adjustment is made to c (i.e. cb).
    xs =real(ifft(fx));
    [si, ciph(iter,:), phi(iter,:)]=cmtmPH(xs,ys,dt,NW);
  end;

  %sorting and averaging to determine confidence levels.
  pl=round(.95*iter);  
  ci=sort(ci);    
  ci=mean(ci(pl,:));
  ci=mean(ci)*ones(size(ci));
  pl=round(.975*iter);  
  phi=sort(phi);  
  phi=[phi(pl,:); -phi(iter-pl+1,:)];
  phi=[mean(phi); -mean(phi)];
  phi(1,1:3)=mean(phi(1,1:3));
  phi(2,1:3)=mean(phi(2,1:3));
  phi(1,end-2:end)=mean(phi(1,end-2:end));
  phi(2,end-2:end)=mean(phi(2,end-2:end));
  temp=conv(phi(1,:),[.25 .5 .25]);
  phi(1,4:end-3)=temp(5:end-4);
  temp=conv(phi(2,:),[.25 .5 .25]);
  phi(2,4:end-3)=temp(5:end-4);
end;

%Cut to one-sided funtions
c = c(pls);
s = s(pls);
ph=ph(pls);

%plotting
if qplot==1,
  %coherence
  figure(gcf); 
  subplot(211); hold on;
  plot(s,c);
  h=ylabel('coherence'); font(h,15);
  if confn>0;
    plot(si,ci*ones(size(si)),'k:');
    pl=find(c>ci(1));
    %title([num2str(100*length(pl)/length(c),2),'% of estimates above 95% confidence level']);
  end;
  axis tight; h=axis; axis([h(1:2) 0 1.025]); h=axis;
  w  = NW/(dt*N);   %half-bandwidth of the dpss
  plot([s(1) h(2)],[1.02 1.02],'k');
  for ds=min(s):2*w:max(s);
    plot([ds ds],[.98 1.02],'k');
  end;
  %axis([h(1) 4 h(3) 1.025]);
  %set(gca,'xscale','log');
  %plot([1/1.5 1/1.5],h(3:4),'k--');
  %plot([1/10 1/10],h(3:4),'k--');
  font(gca,12);
  %h=text(1/70,.15,'(a)'); font(h,15);
  
  %phase
  subplot(212); hold on;
  plot(s,ph);
  if confn>0,
    cu=ph'+phi(1,:)';
    cl=ph'+phi(2,:)';
    c=[.9 .9 .9];
    h=fill([s(1); s(1:end); flipud([s(1:end); s(end)])],[cu(1); cl; flipud([cu; cl(end)])],c);
    set(h,'edgecolor',c);

    pl=find(cu<=180); cu(pl)=-180;
    pl=find(cu> 180); cu(pl)=cu(pl)-360;
    clt=-180*ones(size(cl));
    h=fill([s(1); s(1:end); flipud([s(1:end); s(end)])],[cu(1); clt; flipud([cu; clt(end)])],c);
    set(h,'edgecolor',c);

    pl=find(cl>=-180); cl(pl)=180;
    pl=find(cl< -180); cl(pl)=cl(pl)+360;
    cut=180*ones(size(cl));
    h=fill([s(1); s(1:end); flipud([s(1:end); s(end)])],[cut(1); cl; flipud([cut; cl(end)])],c);    
    set(h,'edgecolor',c);
        
    h=plot(s,ph); %set(h,'linewidth',1.5);      
  end;
  pl=find(s<1/1);
  %plot(s(pl),90*ones(size(pl))+angle(exp(2*i*pi*s(pl)*.21))*180/pi,'k:');
  plot(s,zeros(size(s)),'k:');
  axis tight; h=axis; axis([h(1:2) -180 180]);
  %plot([1/1.5 1/1.5],h(3:4),'k--');
  %plot([1/10 1/10],h(3:4),'k--');
  %set(gca,'xscale','log');
  %xlabel('frequency (cycles/deltat)')
  h=ylabel('phase'); font(h,15);
  font(gca,12);
  %h=text(1/70,-140,'(b)'); font(h,15);
end;









