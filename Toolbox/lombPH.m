%Pn is not scaled by the variance here
%[Power, Prob, freq] = lomb(time, y, freq, linetype)
%
%Lomb's method is used to compute a periodogram where the sum(Power) = var(y).
%
%      freq = frequencies that Power will be computed for; 
%             if length(freq) == 1 and is neg. it is interpreted as delta t 
%             and a freq axis as in fftPH is generated
%      time = times of the measurements y.
%      y    = observations
%      Prob = the probability that the null hypotheses is valid.
%  linetype = what do you want the plot to look like, [] for no plot

function [Pn, Prob, f] = lombPH(t, y, f, linetype)


%		check inputs
deltat = mean(diff(t));
if (length(f)==1 & f(1) < 0), 
df = 1/(length(y)*deltat);
f = df:df:df*fix(length(y)/2);
end;

if length(t) ~= length(y); error('t and y not same length');
	exit; end;

%	detrend y, initialize Pn
%z = detrend(y);
z=y;
N=length(f);
Pn=zeros(size(f));

%	now do main loop for all frequencies
for i=1:length(f)
    w=2*pi*f(i);
    if w > 0 
       twt = 2*w*t;
       tau = atan2(sum(sin(twt)),sum(cos(twt)))/2/w;
       wtmt = w*(t - tau);
       Pn(i) = (sum(z.*cos(wtmt)).^2)/sum(cos(wtmt).^2) + ...
		(sum(z.*sin(wtmt)).^2)/sum(sin(wtmt).^2);
     else
	Pn(i) = (sum(z.*t).^2)/sum(t.^2);
     end
end

%	and normalize by variance, compute probs
PnN=Pn/var(y);
Prob = 1-(1-exp(-PnN)).^N;
for i=1:length(Pn)		% accomodate possible roundoff error
    if Prob(i) < .001
	Prob(i) = N*exp(-PnN(i));
    end
   end
%Pn= Pn*var(y)/sum(Pn);


if length(linetype) > 0  
   plot(f,Pn,linetype)
   title('Lomb Periodogram');
   xlabel(['Cycles/units']);
   ylabel('Power Density, units^2/cycle/deltat');
   set(gca,'yscale','log');   

   if linetype == 2,
   hold on
   place = find(Prob<=.05);   
   plot(f(place),Pn(place),'*')

   for count = 1:length(place)
   h=text(f(place(count)),Pn(place(count)),[' ',num2str(1/f(place(count)),3)]);
   set(h,'fontsize',8);
   end;
   end;      
axis tight
end;

%display('Lomb: Sum of Power, Variance of Ydata');
%display([sum(Pn), var(y)]);
