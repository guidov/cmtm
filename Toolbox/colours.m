% colours.m
% defines a bunch of colour maps 
% written by Camille Li, UW 2003.

% thermal1 (blue white yellow)
%-----------------------------
% this is my version of dave thompson's colour map
% various other thermals follow - they're all similar
% thermalline is good for contour plots
cmap = [ ...
10 50 120;
15 75 165;
30 110 200;
60 160 240;
80 180 250;
130 210 255;
160 240 255;
200 250 255;
230 255 255;
255 255 255;
255 255 255;
255 250 220;
255 232 120;
255 192 60;
255 160 0;
255 96 0;
255 50 0;
225 20 0;
192 0 0;
165 0 0];
thermal1 = (cmap./255);

r1=[0.55 0.2  0.5 0.8  1.0 1.0 1.0  1.0 1.0 0.78];
g1=[0.25 0.5  0.8 1.0  1.0 1.0 1.0  0.7 0.5 0.15];
b1=[0.75 0.84 1.0 0.95 1.0 1.0 0.25 0   0   0];
thermal2=[r1',g1',b1'];

thermal3=thermal2(2:end,:);
thermal3=[[0.2 0 0.2];
          [0 0 0.35];
          [0.1 0.25 0.5];
          thermal3;
          [0.5 0 0];
          [0.3 0 0]];


	
	
thermal4=[thermal1(1:7,:); [1 1 1]; thermal1(14:20,:)];
thermal5=[thermal1(1:7,:); thermal1(14:20,:)];





return;

thermalline=[0.0392    0.1961    0.4706;
             0.3137    0.7059    0.9804;
             0.6275    0.9412    1.0000;
             0.6275    1.0000    0.2553;
             1.0000    0.7529    0.2353;
	     1.0000    0.3765	 0;			                 
	     0.6471    0	 0];


colormap(thermal4); 
clear thermal* cmap thermalline;


% rwb (red white blue)
%---------------------
cmin=-1; cmax=1;
r1=[((0:31))/31,ones(1,31)];
g1=[((0:31))/31,((30:-1:0))/31];
b1=[ones(1,31),((31:-1:0))/31];
rwb=[r1',g1',b1'];

% algae (blue white green)
%-------------------------
r1=[zeros(1,21),(0:0.1:1),(.9:-.1:0),zeros(1,21)];
g1=[((1:31))/31,linspace(1,0.5,32)];
b1=[ones(1,31),(linspace(31,0.1,32))/31];
algae=[r1',g1',b1'];

% burnt (like rwb but "smokier")
%-------------------------------
r1=[linspace(0.3,1,25),ones(1,6),...
    1,[0.9 0.8 0.7 0.6],linspace(0.5,0,27)];
g1=[linspace(0,0.7,27),[0.75 0.8 0.85 0.9],...
    1,[0.9 0.85 0.8 0.7],linspace(0.7,0,27)];
b1=[linspace(0.1,0.5,27),[0.6 0.7 0.8 0.9]...
    1,[1 1 1 1],linspace(1,0.3,27)];
burnt=[r1',g1',b1'];
