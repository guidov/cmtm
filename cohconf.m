%function [cl]=cohconf(v,level,unbias,c);
%
%inputs:   v     - degrees of freedom
%          level - confidence level, i.e. 0.95 gives the 95% c.l. (.95=default).
%          bias  - are coherence estimates bias corrected (use cohbias.m), 0=no (default), 1=yes.
%          c     - true coherence, 0 <= c < 1 (0=default).
%                    
%
%outputs:  
%          cl     - coherence value for selected confidence level 
%
%
%
%Peter Huybers
%phuybers@mit.edu
%MIT, 2003

function [cl]=cohconf(n,level,unbias,c);

if nargin<1, help cohconf; return; end;
if nargin<2, level=.95; end;
if nargin<3, unbias=0;    end;
if nargin<4, c=0;       end;

if n<2, disp('Warning: degress of freedom must be greater or equal to two.'); return;   end;
if level<=0 | level>=1, disp('Warning: confidence level should be between zero and one.'); return; end;

%Calculated according to: Amos and Koopmans, "Tables of the distribution of the
%coefficient of coherence for stationary bivariate Gaussian processes", Sandia
%Corporation, 1963
%Also see Priestly, 1981
z=0:.0005:1;
for i1=1:length(z),
  A(1)=1;
  for k=1:n-1;
    A(k+1)=A(k)*(n-k)*(2*k-1)/((2*n-(2*k+1))*k)*((1-c*z(i1))/(1+c*z(i1)))^2;
  end;
  f(i1)=2*(n-1)*(1-c^2)^n*z(i1)*(1-z(i1)^2)^(n-2)/((1+c*z(i1))*(1-c*z(i1))^(2*n-1))...
	*gamma(n-.5)/(sqrt(pi)*gamma(n))*sum(A);
end;


%Use a quadratic Newton-Cotes methods to determine the cumulative sum
for i1 = 2:length(f)/2;
  F(i1) = f(2*(i1-1)+1) + 4*f(2*i1) + f(2*i1+1);
end

F = F/(6*length(F));
F=[fliplr(1-cumsum(fliplr(F))/sum(F)) 1];
Fz=[z(1:2:end-2) 1]; 
pl=find(diff(F)>0); pl=[1 pl+1];
cl=interp1(F(pl),Fz(pl),level);

if unbias==1,
  cl=cohbias(n,cl);
end;
