%A matrix smoother. 
%
%function [SM] = smooth2PH(M,dim,iter);
%
%SM   = smoothed output Matrix.
%M    = input Matrix.
%dim  = dimension of smoother, must be an odd number
%iter = self-convolutions of the smoother, 0=boxcar.  

function [SM] = smooth2PH(M,dim,iter);

if ~exist('dim'),  dim=3;  end;
if ~exist('iter'), iter=1;  end;
if rem(dim,2)~=1, disp('??? Error using smooth2PH ==> smoother dimension must be odd'); SM=NaN; return; end;

x=ones(dim,dim);
for ct=1:iter,
  x=conv2(x,x); 
end;
x=x/sum(sum(x));
SM=conv2(x,M)./conv2(x,ones(size(M)));

n=(length(x)-1)/2;
SM=SM(n+1:end-n,n+1:end-n);

