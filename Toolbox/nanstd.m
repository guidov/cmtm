function sx = nanstd(x)

%Returns standard deviation excluding the nans.  

  % if x is a vector, make sure it is a row vector
  if length(x)==prod(size(x))         
    x = x(:);                         
  end  
  [m,n]   = size(x);

  for ct=1:n;
    pl = find(isnan(x(:,ct))==0);
    sx(ct) = std(x(pl,ct)); 
  end;
 
