%function  [s, c, ph, ci, phi] = cmtm(x,y,dt,NW,qbias,confn,qplot);
%
%Multi-taper method coherence using adaptive weighting and correcting 
%for the bias inherent to coherence estimates.  The 95% coherence
%confidence level is computed by cohconf.m.  In addition, a built-in
%Monte Carlo estimation procedure is available to estimate phase 95%
%confidence limits. 
%
% Inputs:
%         x     - Input data vector 1.  
%         y     - Input data vector 2.  
%         dt    - Sampling interval (default 1) 
%         NW    - Number of windows to use (default 8) 
%         qbias - Correct coherence estimate for bias (yes, 1)  (no, 0, default).
%         confn - Number of iterations to use in estimating phase uncertainty using a Monte Carlo
%                 method. (default 0)
%         qplot - Plot the results, (yes, 1), (No, 0, default).  The upper tickmarks indicate the
%                 bandwidth of the coherence and phase estimates.  
%
% Outputs:
%         s       - frequency
%         c       - coherence
%         ph      - phase
%         ci      - 95% coherence confidence level
%         phi     - 95% phase confidence interval, bias corrected
%                   (add and subtract phi from ph).
%
%
%required files: cohconf.m, cohbias.m, cohbias.mat, Matlab signal processing toolbox.
%
%Peter Huybers
%MIT, 2003
%phuyber@mit.edu

function [s, c, ph, ci, phi] = cmtm(x,y,dt,NW,qbias,confn,qplot);

%check input
if nargin<2, help cmtm; return; end;
if nargin<7, qplot=0;  end;
if nargin<6, confn=0;  end;
if nargin<5, qbias=0;  end;
if nargin<4, NW=8;     end;
if length(NW)==0, NW=8;end;
if nargin<3, dt=1;     end;
if length(dt)==0, dt=1;end;


if NW<1.5, disp('Warning: NW must be greater or equal to 1.5'); return; end;
if nargin>4, 
	  disp('-------------------------');
disp(['Number of windows: ',num2str(NW)]);
if qbias==1, disp('Bias correction:   On');
else,        disp('Bias correction:   Off'); end;
disp(['Confidence Itera.: ',num2str(confn)]);
if qplot==1, disp('Plotting:          On');
else,        disp('Plotting:          Off'); end;
disp('-------------------------');
end;

x=x(:)-mean(x); 
y=y(:)-mean(y);
if length(x)~=length(y), disp('Warning: the lengths of x and y must be equal.'); return; end;

%define some parameters
N   = length(x);
k   = min(round(2*NW),N); 
k   = max(k-1,1);
s   = (0:1/(N*dt):1/dt-1/(N*dt))';
pls=2:(N+1)/2+1;
v   = (2*NW-1); %approximate degrees of freedom

if rem(length(y),2)==1; pls=pls(1:end-1); end;

%Compute the discrete prolate spheroidal sequences, requires the spectral analysis toolbox.
[E,V]=dpss(N,NW,k);

%Compute the windowed DFTs.
fkx=fft(E(:,1:k).*x(:,ones(1,k)),N);
fky=fft(E(:,1:k).*y(:,ones(1,k)),N);

Pkx=abs(fkx).^2;
Pky=abs(fky).^2;

%Itteration to determine adaptive weights:    
for i1=1:2,
  if i1==1,   vari=x'*x/N; Pk=Pkx; end;
  if i1==2,   vari=y'*y/N; Pk=Pky; end;
  P    = (Pk(:,1)+Pk(:,2))/2;   % initial spectrum estimate
  Ptemp= zeros(N,1);
  P1   = zeros(N,1);
  tol  = .0005*vari/N;          % usually within 'tol'erance in about three iterations, see equations from [2] (P&W pp 368-370).   
  a    = vari*(1-V);
  while sum(abs(P-P1)/N)>tol            
    b=(P*ones(1,k))./(P*V'+ones(N,1)*a'); % weights
    wk=(b.^2).*(ones(N,1)*V');            % new spectral estimate
    P1=(sum(wk'.*Pk')./ sum(wk'))';
    Ptemp=P1; P1=P; P=Ptemp;              % swap P and P1
  end
  if i1==1, 
    fkx=sqrt(k)*sqrt(wk).*fkx./repmat(sum(sqrt(wk'))',1,k);
    Fx=P;  %Power density spectral estimate of x
  end;
  if i1==2, 
    fky=sqrt(k)*sqrt(wk).*fky./repmat(sum(sqrt(wk'))',1,k);
    Fy=P;  %Power density spectral estimate of y
  end;
end;
%As a check, the quantity sum(abs(fkx(pls,:))'.^2) is the same as Fx and
%the spectral estimate from pmtmPH.

%Compute coherence
Cxy= sum([fkx.*conj(fky)]');
ph = angle(Cxy)*180/pi;
c  = abs(Cxy)./sqrt(sum(abs(fkx').^2).*sum(abs(fky').^2));

%correct for the bias of the estimate
if qbias==1,
  c=cohbias(v,c)'; 
end;

%Phase uncertainty estimates via Monte Carlo analysis. 
if confn>1,
  cb=cohbias(v,c)';
for iter=1:confn;
  if rem(iter,10)==0, disp(['phase confidence iteration: ',num2str(iter)]); end;
  fx=fft(randn(size(x))+1);
  fx=fx/sum(abs(fx));
  fy=fft(randn(size(y))+1); 
  fy=fy/sum(abs(fy));
  ys =real(ifft(fy.*sqrt(1-cb'.^2)));    
  ys =ys+real(ifft(fx.*cb')); 
  xs =real(ifft(fx));
  [si, ciph(iter,:), phi(iter,:)]=cmtm(xs,ys,dt,NW);
end;
  pl=round(.975*iter);     %sorting and averaging to determine confidence levels.
  phi=sort(phi);  
  phi=[phi(pl,:); -phi(iter-pl+1,:)];
  phi=mean(phi);
  phi=conv(phi(1:end),[1 1 1]/3);  
  phi=phi(2:end-1);
else,
  phi=zeros(size(pls)); 
end;
  
%Cut to one-sided funtions
c = c(pls);
s = s(pls)';
ph=ph(pls);
phl=ph-phi;
phu=ph+phi;

%Coherence confidence level
ci=cohconf(v,.95);  %not corrected for bias, this is conservative.
ci=ci*ones(size(c));

%plotting
if qplot==1,
  %coherence
  figure(gcf); clf;
  subplot(211); hold on;
  plot(s,c);
  h=ylabel('coherence');
  h=xlabel('frequency');
  plot(s,ci,'k--');
  pl=find(c>ci(1));
  title(['mean is ',num2str(mean(c),2),'   ',num2str(100*length(pl)/length(c),2),'% of estimates above 95% confidence level']);
  axis tight; h=axis; axis([h(1:2) 0 1.025]);
  w  = NW/(dt*N);   %half-bandwidth of the dpss
  plot([s(1) h(2)],[1.02 1.02],'k');
  for ds=min(s):2*w:max(s);
    plot([ds ds],[.98 1.02],'k');
  end;
  
  %phase
  subplot(212); hold on;
  plot(s,ph);
  if confn>0,
    col=[.9 .9 .9];
    h=fill([s(1) s(1:end) fliplr([s(1:end) s(end)])],[phu(1) phl fliplr([phu phl(end)])],col);
    set(h,'edgecolor',col);

    pl=find(phu<=180); phu(pl)=-180;
    pl=find(phu> 180); phu(pl)=phu(pl)-360;
    phlt=-180*ones(size(phl));
    h=fill([s(1) s(1:end) fliplr([s(1:end) s(end)])],[phu(1) phlt fliplr([phu phlt(end)])],col);
    set(h,'edgecolor',col);

    pl=find(phl>=-180); phl(pl)=180;
    pl=find(phl< -180); phl(pl)=phl(pl)+360;
    phut=180*ones(size(phl));
    h=fill([s(1) s(1:end) fliplr([s(1:end) s(end)])],[phut(1) phl fliplr([phut phl(end)])],col);    
    set(h,'edgecolor',col);        
  end;
  h=plot(s,ph); 
  plot(s,zeros(size(s)),'k--');
  axis tight; h=axis; axis([h(1:2) -180 180]);
  h=xlabel('frequency'); 
  h=ylabel('phase'); 
end;



