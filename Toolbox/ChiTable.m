%This makes a tabel of Chi squared distributions.

clear; 
v=[1:5]; Per=[.55:.005:.99];  ct1=0;
for c1=Per;
  ct2=0;
  ct1=ct1+1;
  c1
  2*c1-1
for c2=v;
  ct2=ct2+1;
chi=chi2conf(2*Per-1,c2*2);
level(ct2,ct1)=chi(2); 

end; end;