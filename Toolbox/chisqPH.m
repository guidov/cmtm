function p = chisqPH(x,p);
%function p = chisqPH(x,p);
%
%Chi squared confidence level based on value of x and p independent variables.
%This is calculated witht the modified gamma function with A = p/2 and B =2
%from code taken within the chi2conf function given in matlab

%GAMCDF Gamma cumulative distribution function.
%   P = GAMCDF(X,A,B) returns the gamma cumulative distribution
%   function with parameters A and B at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   Some references refer to the gamma distribution with a single
%   parameter. This corresponds to the default of B = 1. 
%
%   GAMMAINC does computational work.

%   References:
%      [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%      Springer-Verlag, 1986. p. 401.
%      [2]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.32.

%   Was: Revision: 1.2, Date: 1996/07/25 16:23:36

a = p/2;
b = 2;

%   Return NaN if the arguments are outside their respective limits.
k1 = find(a <= 0 | b <= 0);     
if any(k1)
    tmp   = NaN;
    p(k1) = tmp(ones(size(k1)));
end

% Initialize P to zero.
p = zeros(size(x));

k = find(x > 0 & ~(a <= 0 | b <= 0));
if any(k), 
    p(k) = gammainc(x(k) ./ b(k),a(k));
end

% Make sure that round-off errors never make P greater than 1.
k = find(p > 1);
if any(k)
    p(k) = ones(size(k));
end


function y = gampdf(x,a,b)
%GAMPDF Gamma probability density function.
%   Y = GAMPDF(X,A,B) returns the gamma probability density function 
%   with parameters A and B, at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   Some references refer to the gamma distribution with a single
%   parameter. This corresponds to the default of B = 1.

%   References:
%      [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%      Springer-Verlag, 1986, pages 401-402.

%   Was: Revision: 1.2, Date: 1996/07/25 16:23:36

if nargin < 3, 
    b = 1; 
end

if nargin < 2, 
    error('Requires at least two input arguments'); 
end

[errorcode x a b] = distchck(3,x,a,b);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

% Initialize Y to zero.
y = zeros(size(x));

%   Return NaN if the arguments are outside their respective limits.
k1 = find(a <= 0 | b <= 0);     
if any(k1)
    tmp = NaN;
    y(k1) = tmp(ones(size(k1)));
end

k=find(x > 0 & ~(a <= 0 | b <= 0));
if any(k)
    y(k) = (a(k) - 1) .* log(x(k)) - (x(k) ./ b(k)) - gammaln(a(k)) - a(k) .* log(b(k));
    y(k) = exp(y(k));
end
k1 = find(x == 0 & a < 1);
if any(k1)
  tmp = Inf;
  y(k1) = tmp(ones(size(k1)));
end
k2 = find(x == 0 & a == 1);
if any(k2)
  y(k2) = (1./b(k2));
end

function [errorcode,out1,out2,out3,out4] = distchck(nparms,arg1,arg2,arg3,arg4)
%DISTCHCK Checks the argument list for the probability functions.

%   B.A. Jones  1-22-93
%   Was: Revision: 1.2, Date: 1996/07/25 16:23:36

errorcode = 0;

if nparms == 1
    out1 = arg1;
    return;
end
    
if nparms == 2
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    scalararg1 = (prod(size(arg1)) == 1);
    scalararg2 = (prod(size(arg2)) == 1);
    if ~scalararg1 & ~scalararg2
        if r1 ~= r2 | c1 ~= c2
            errorcode = 1;
            return;         
        end     
    end
    if scalararg1
        out1 = arg1(ones(r2,1),ones(c2,1));
    else
        out1 = arg1;
    end
    if scalararg2
        out2 = arg2(ones(r1,1),ones(c1,1));
    else
        out2 = arg2;
    end
end
    
if nparms == 3
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    [r3 c3] = size(arg3);
    scalararg1 = (prod(size(arg1)) == 1);
    scalararg2 = (prod(size(arg2)) == 1);
    scalararg3 = (prod(size(arg3)) == 1);

    if ~scalararg1 & ~scalararg2
        if r1 ~= r2 | c1 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg1 & ~scalararg3
        if r1 ~= r3 | c1 ~= c3
            errorcode = 1;
            return;                 
        end
    end

    if ~scalararg3 & ~scalararg2
        if r3 ~= r2 | c3 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg1
      out1 = arg1;
    end
    if ~scalararg2
      out2 = arg2;
    end
    if ~scalararg3
      out3 = arg3;
    end
    rows = max([r1 r2 r3]);
   columns = max([c1 c2 c3]);  
       
    if scalararg1
        out1 = arg1(ones(rows,1),ones(columns,1));
    end
   if scalararg2
        out2 = arg2(ones(rows,1),ones(columns,1));
   end
   if scalararg3
       out3 = arg3(ones(rows,1),ones(columns,1));
   end
     out4 =[];
    
end

if nparms == 4
    [r1 c1] = size(arg1);
    [r2 c2] = size(arg2);
    [r3 c3] = size(arg3);
    [r4 c4] = size(arg4);
    scalararg1 = (prod(size(arg1)) == 1);
    scalararg2 = (prod(size(arg2)) == 1);
    scalararg3 = (prod(size(arg3)) == 1);
    scalararg4 = (prod(size(arg4)) == 1);

    if ~scalararg1 & ~scalararg2
        if r1 ~= r2 | c1 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg1 & ~scalararg3
        if r1 ~= r3 | c1 ~= c3
            errorcode = 1;
            return;                 
        end
    end

    if ~scalararg1 & ~scalararg4
        if r1 ~= r4 | c1 ~= c4
            errorcode = 1;
            return;                 
        end
    end

    if ~scalararg3 & ~scalararg2
        if r3 ~= r2 | c3 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg4 & ~scalararg2
        if r4 ~= r2 | c4 ~= c2
            errorcode = 1;
            return;         
        end
    end

    if ~scalararg3 & ~scalararg4
        if r3 ~= r4 | c3 ~= c4
            errorcode = 1;
            return;         
        end
    end


    if ~scalararg1
       out1 = arg1;
    end
    if ~scalararg2
       out2 = arg2;
    end
    if ~scalararg3
      out3 = arg3;
    end
    if ~scalararg4
      out4 = arg4;
    end
 
   rows = max([r1 r2 r3 r4]);
   columns = max([c1 c2 c3 c4]);     
    if scalararg1
       out1 = arg1(ones(rows,1),ones(columns,1));
   end
   if scalararg2
       out2 = arg2(ones(rows,1),ones(columns,1));
   end
   if scalararg3
       out3 = arg3(ones(rows,1),ones(columns,1));
   end
   if scalararg4
       out4 = arg4(ones(rows,1),ones(columns,1));
   end
end

