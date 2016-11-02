%Find half-height of the terminations
%
%function [it]=HalfHeigth(t,y,it);

function [it]=HalfHeight(t,y,it);

t=t(:);
y=y(:);

%y=smoothPH(y,5);
for ct=1:length(it); 
  [dum pl]=min(abs(t-it(ct)));
  plm=pl;
  
  while (y(plm)>y(plm+1) & plm>2)
     plm=plm-1; 
     pl=plm;
  end;
  
  
  while (y(plm-1)<y(plm) & plm>2),
    plm=plm-1;
  end;
  plp=pl;
  while (y(plp+1)>y(plp) & plp<length(y)-1),
    plp=plp+1;
  end;
  hh=(y(plp)+y(plm))/2;
  if plm<pl & plp>pl, 
    it(ct)=interp1(y(plm:plp),t(plm:plp),hh);
  end;
end;
