function t=chiaux(chisq,vec)
% Aux. function used by CHISQUARED_TABLE
%
%Peter R. Shaw, WHOI

% Peter R. Shaw, Woods Hole Oceanographic Institution
% Woods Hole, MA 02543   pshaw@whoi.edu
% March, 1990;  
% Last Revision: Oct 1992 gamma/incgamma with version 4

nu=vec(1);
P=vec(2);
versn_str=version; eval(['versn=' versn_str(1) ';']);
if versn<=3,
  t = P-gamma(nu/2,chisq/2);
else
  t = P-gammainc(chisq/2,nu/2);
end
