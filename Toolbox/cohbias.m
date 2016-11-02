%function [cu]=cohbias(v,cb);
%
%Corrects for the bias inherent to coherence estimates.  Note the Matlab
%function cohere.m returns squared-coherence, and the square-root should
%be used.  Coherence below the minimum expected value returns a zero. 
%
%Requires the file cohbias.mat.  If the file does not exist,
%prompts whether it should be created -- note the calculation
%takes roughly an hour on a 2 GHz machine (i.e. it should be 
%easier to get the file from http://web.mit.edu/~phuybers/www/XCM/index.html.)    
%
%
%inputs:   v   - degrees of freedom, single value or vector (2 <= n <= 50)
%          cb - biased coherence estimate, single value of vector (0 <= c <= 1).
%          
%outputs:  cu  - unbiased cohernce estimate (always less than cb).
%
%
%Peter Huybers
%phuybes@mit.edu
%MIT, 2003.

function [cu]=cohbias(v,cb);

if nargin<2, help cohbias; return; end;

if v<2, disp('Warning: degress of freedom must be greater or equal to two.'); return;   end;
if cb<0 | cb>1, disp('Warning: biased coherence should be between zero and one, inclusive.'); return; end;


if v>50, disp('using 50 degrees of freedom'); v=50; end;
if nargin==0; help cohbias.m; return; end;

if exist('cohbias.mat')==0;
  %Cohbias.mat file should be down-loaded with cmtm.m and cohbias.m
  %The routine is included primarily to show how it was created.
  n=2:1:50;
  c=.1:.0001:1;
  disp('-- The file cohbias.mat does not exist within the path.');
  qans=input('-- To create this file now enter ''y'' or to skip ''n''. \n--> ','s');
  qans,
  if strncmpi(qans,'y',1);	      
    z=0:.1:1;
    for i3=1:length(n);
      disp(n(i3)),
      for i2=1:length(c)-1;
	for i1=1:length(z),
	  A(1)=1;
	  %Calculated according to: Amos and Koopmans, "Tables of the distribution of the
	  %coefficient of coherence for stationary bivariate Gaussian processes", Sandia
	  %Corporation, 1963
	  %
	  %Also see the manuscript of Wunsch, C. "Time-Series Analysis.  A Heuristic Primer".  
	  for k=1:n(i3)-1;
	    A(k+1)=A(k)*(n(i3)-k)*(2*k-1)/((2*n(i3)-(2*k+1))*k)*((1-c(i2)*z(i1))/(1+c(i2)*z(i1)))^2;
	  end;
	  f(i1)=2*(n(i3)-1)*(1-c(i2)^2)^n(i3)*z(i1)*(1-z(i1)^2)^(n(i3)-2)/((1+c(i2)*z(i1))*(1-c(i2)*z(i1))^(2*n(i3)-1))...
		*gamma(n(i3)-.5)/(sqrt(pi)*gamma(n(i3)))*sum(A);
	end;
	%Use a quadratic Newton-Cotes methods to determine the cumulative sum
	for i1 = 2:length(f)/2;
	  M(i1) = [f(2*(i1-1)+1) + 4*f(2*i1) + f(2*i1+1)]*z(2*i1);
	end
	expect(i3,i2)=sum([M 1])/(6*(length(M)));
      end;
      expect(i3,i2+1)=1;
    end;
    save cohbias.mat expect n c; 
  else, %if skip cohbias.mat calculation
    expect=repmat(c,length(n),1);
  end;  %stop qans condition
else,   %if cohbias.mat already exists
  load cohbias.mat expect n c;
end;    %stop cohbias calculation

cb=cb(:);
c=c(:);
n=n(:);
v=v(:);

for ct=1:length(c);
  ec(ct)=interp1(n,expect(:,ct),v);
end;

for ct=1:length(cb);
  cu(ct)=interp1(ec,c,cb(ct));
end;

cu=cu(:);

pl=find(isnan(cu)==1 & cb<1 & cb>=0); %If cu is NaN while cb is between (0,1)
cu(pl)=0; 




