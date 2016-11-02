function h=circle(r,x0,y0,C,Nb)
% CIRCLE adds circles to the current plot
%
% CIRCLE(r,x0,y0) adds a circle of radius r centered at point x0,y0. 
% If r is a vector of length L and x0,y0 scalars, L circles with radii r 
% are added at point x0,y0.
% If r is a scalar and x0,y0 vectors of length M, M circles are with the same 
% radius r are added at the points x0,y0.
% If r, x0,y0 are vector of the same length L=M, M circles are added. (At each
% point one circle).
% if r is a vector of length L and x0,y0 vectors of length M~=L, L*M circles are
% added, at each point x0,y0, L circles of radius r.
%
% CIRCLE(r,x0,y0,C)
% adds circles of color C. C may be a string ('r','b',...) or the RGB value. 
% If no color is specified, it makes automatic use of the colors specified by 
% the axes ColorOrder property. For several circles C may be a vector.
%
% CIRCLE(r,x0,y0,C,Nb), Nb specifies the number of points used to draw the 
% circle. The default value is 300. Nb may be used for each circle individually.
%
% h=CIRCLE(...) returns the handles to the circles.
%
% Try out the following (nice) examples:
%
%% Example 1
%
% clf;
% x=zeros(1,200);
% y=cos(linspace(0,1,200)*4*pi);
% rad=linspace(1,0,200);
% cmap=hot(50);
% circle(rad,x,y,[flipud(cmap);cmap]);
%
%% Example 2
%
% clf;
% the=linspace(0,pi,200); 
% r=cos(5*the);
% circle(0.1,r.*sin(the),r.*cos(the),hsv(40));
% 
%
%% Example 3
%
% clf
% [x,y]=meshdom(1:10,1:10);
% circle([0.5,0.3,0.1],x,y,['r';'y']);
%
%% Example 4
%
% clf
% circle(1:10,0,0,[],3:12);
%
%% Example 5
%
% clf;
% circle((1:10),[0,0,20,20],[0,20,20,0]);

% written by Peter Blattner, Institute of Microtechnology, University of 
% Neuchatel, Switzerland, blattner@imt.unine.ch



% Check the number of input arguments 

if nargin<1,
  r=[];
end;

if nargin==2,
  error('Not enough arguments');
end;

if nargin<3,
  x0=[];
  y0=[];
end;
 

if nargin<4,
  C=[];
end

if nargin<5,
  Nb=[];
end

% set up the default values

if isempty(r),r=1;end;
if isempty(x0),x0=0;end;
if isempty(y0),y0=0;end;
if isempty(Nb),Nb=300;end;
if isempty(C),C=get(gca,'colororder');end;

% work on the variable sizes

x0=x0(:);
y0=y0(:);
r=r(:);
Nb=Nb(:);


if isstr(C),C=C(:);end;

if length(x0)~=length(y0),
  error('length(x0)~=length(y0)');
end;

% how many rings are plottet

if length(r)~=length(x0)
  maxk=length(r)*length(x0);
else
  maxk=length(r);
end;

% drawing loop

for k=1:maxk
  
  if length(x0)==1
    xpos=x0;
    ypos=y0;
    rad=r(k);
  elseif length(r)==1
    xpos=x0(k);
    ypos=y0(k);
    rad=r;
  elseif length(x0)==length(r)
    xpos=x0(k);
    ypos=y0(k);
    rad=r(k);
  else
    rad=r(fix((k-1)/size(x0,1))+1);
    xpos=x0(rem(k-1,size(x0,1))+1);
    ypos=y0(rem(k-1,size(y0,1))+1);
  end;

  the=linspace(0,2*pi,Nb(rem(k-1,size(Nb,1))+1,:)+1);
  h(k)=line(rad*cos(the)+xpos,rad*sin(the)+ypos);
  set(h(k),'color',C(rem(k-1,size(C,1))+1,:));

end;