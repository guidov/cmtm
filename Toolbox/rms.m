
function [x]=rms(y);

y=y(:);

x=sqrt(nanmean(y.^2));
