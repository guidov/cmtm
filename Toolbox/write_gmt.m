function write_gmt(x,y,z,cdf)
% WRITE_GMT
% writes Matlab array z into GMT .grd file
% 
% Usage:  write_gmt(x,y,z,cdf_file)
% Example:  contour(x,y,z);
%           write_gmt(x,y,z,grdfile);
%           [x,y,z]=read_gmt(grdfile);
%           contour(x,y,z);
%
%  see also READ_GMT

% Rich Signell
% rsignell@usgs.gov
% M. Tivey mod mexcdf60 for MATLAB V6 Oct 31 2001
% BL. Owems mod mexcdf60 for MATLAB V6 Sep 2003
%
z=flipud(z)';
[nx,ny]=size(z);
xysize=nx*ny;
x1=min(x);
x2=max(x);
y1=min(y);
y2=max(y);
spacing=[(x2-x1)/(nx-1) (y2-y1)/(ny-1)];
x_range=[x1 x2];
y_range=[y1 y2];
z_range=[min(z(:)) max(z(:))];
dimension=[nx ny];
%
% Create netCDF file
cdfid = mexcdf53('create',cdf,'CLOBBER');
% Suppress error messages from netCDF
[rcode]= mexcdf53('setopts',0);
%
% Define dimensions
dim1 = mexcdf53('dimdef',cdfid,'side',2);
dim2 = mexcdf53('dimdef',cdfid,'xysize',xysize);
% Define variables
%
x_range_id   = mexcdf53('vardef',cdfid,'x_range','double',1,dim1);
  mexcdf53('attput',cdfid,x_range_id,'units','CHAR',-1,'user_x_unit');
y_range_id   = mexcdf53('vardef',cdfid,'y_range','double',1,dim1);
  mexcdf53('attput',cdfid,y_range_id,'units','CHAR',-1,'user_y_unit');
z_range_id   = mexcdf53('vardef',cdfid,'z_range','double',1,dim1);
  mexcdf53('attput',cdfid,z_range_id,'units','CHAR',-1,'user_z_unit');
spacing_id   = mexcdf53('vardef',cdfid,'spacing','double',1,dim1);
dimension_id   = mexcdf53('vardef',cdfid,'dimension','long',1,dim1);
z_id   = mexcdf53('vardef',cdfid,'z','float',1,dim2);
  mexcdf53('attput',cdfid,z_id,'long_name','CHAR',-1,'surface');
  mexcdf53('attput',cdfid,z_id,'scale_factor','float',1,1.);
  mexcdf53('attput',cdfid,z_id,'add_offset','float',1,0.);
  mexcdf53('attput',cdfid,z_id,'node_offset','long',1,0);
% Create Global attributes
mexcdf53('attput',cdfid,'GLOBAL','title','CHAR',-1,'surface');
mexcdf53('attput',cdfid,'GLOBAL','source','CHAR',-1,'write_gmt (matlab)');
%
% End define mode
mexcdf53('endef',cdfid);

% Store coordinate (independent) variables
mexcdf53('varput',cdfid,x_range_id,0,2,x_range);
mexcdf53('varput',cdfid,y_range_id,0,2,y_range);
mexcdf53('varput',cdfid,z_range_id,0,2,z_range);
mexcdf53('varput',cdfid,spacing_id,0,2,spacing);
mexcdf53('varput',cdfid,dimension_id,0,2,dimension);
n=length(z(:));
mexcdf53('varput',cdfid,z_id,0,n,z(:));
%close netCDF file
mexcdf53('close',cdfid);
disp([cdf ' created'])

