function [x,y,z]=read_gmt(grdfile)
% READ_GMT 
% reads GMT grid file into Matlab array z
%   Usage: [x,y,z]=read_gmt(grdfile);
%    where 
%          x = east coordinate vector (eg. longitude)
%          y = north coordinate vector (eg. latitude)
%          z = matrix of gridded values (eg. bathy grid)
%
%   Example:
%           [x,y,z]=read_gmt('foo.grd');
%           contour(x,y,z)

% Rich Signell
% rsignell@usgs.gov
% M. Tivey mod to MATLAB5 using mexcdf53  May 1997
% M. Tivey mod to MATLAB6 using mexcdf60 Oct 2001
% BL. Owens mod to MATLAB6 using mexcdf53 Sep 2003

cdfid=mexcdf53('open',grdfile,'nowrite');
oldopts=mexcdf53('setopts',0);
x_range=mexcdf53('varget',cdfid,'x_range',0,2);
y_range=mexcdf53('varget',cdfid,'y_range',0,2);
spacing=mexcdf53('varget',cdfid,'spacing',0,2);
dims=mexcdf53('varget',cdfid,'dimension',0,2);
nx=dims(1);
ny=dims(2);
xysize=nx*ny;
z=mexcdf53('varget',cdfid,'z',0,xysize);
mexcdf53('close',cdfid);
z=reshape(z,nx,ny);
z=flipud(z.');
x=x_range(1):spacing(1):x_range(2);
y=y_range(1):spacing(2):y_range(2);

