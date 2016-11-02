

%This function makes a given matrix E into a square symmetric matrix 
% of the form [0 E'; E 0]


function [sqE] = sqm(E);

vl = length(E(:,1));
hl = length(E(1,:));
sqE  = zeros(vl+hl);
sqE(1:hl,(hl+1):(hl+vl)) = E';
sqE(hl+1:hl+vl,1:hl)     = E;
