
%This function takes data with times specified at certain depths and makes linear time line
%between those specified times.  It also interpolates the times to step kyrs.  
%
%        function [ot,od] = mkb1(it, id, depth, idata, step) 
%
% ot = output time scale, evenly space at step kyrs
% od = output data that corresponds to that time
% it = input time marks
% id = corresponding input depth marks
% depth = the depth that each measurement in idata was taken at
% idata = the data we are developing the new time scale for
% step  = the kyrs between each ot.  

function [ot,od] = interpPH(it, id, depth, idata, step) 

if ~exist('it') | ~exist('id') | ~exist('depth') | ~exist('idata') | ~exist('step')
   disp('Some input parameters are missing; using default values')
   it = [0 780 900 1020 1700 1870 2430];
   id = [0 3180 4026 4317 7355 8141 11130];
   m607
   depth = C607(:,1);
   idata = C607(:,2);
   clear m607
   step = 4;   
end;

%otrough(1:length(depth)) = [];
for place = 1:(length(id)-1)
  t1 = it(place);
  t2 = it(place+1);
  d1 = id(place);
  d2 = id(place+1);
  [d,pd1] =min(abs(depth(:,1)-d1));
  [d,pd2] =min(abs(depth(:,1)-d2));
  if place+1 == length(id), pd2 = length(depth); end;
  %Linear Interpolation: (percent of depth) * (time interval) + (time begin)
  for row = pd1:pd2
    otrough(row) =  (depth(row,1)-d1)/(d2-d1) * (t2-t1) + t1; 
  end; %stop interpolating loop
end; %stop magnetic time/depth loop

%Making time steps evenly spaced
ot = 0:step:otrough(length(otrough))';
od = interp1(otrough,idata,ot);
