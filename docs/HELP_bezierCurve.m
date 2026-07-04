%% bezierCurve
% Below is a demonstration of the features of the |bezierCurve| function

%%
clear; close all; clc;

%% Syntax
% |[V]=bezierCurve(p,n);|

%% Description
% This function uses the control points p to create a Bézier curve. The
% output V consists of n points on the curve. 

%% Examples

%%
% Plot settings

fontSize=20;
lineWidth1=2;
lineWidth2=4;
markerSize1=60;
markerSize2=40;

%% Example 1: 

p=[0 0 0; 1 1 1; 2 -1 1; 3 0 0]; %Control points
n=25; %Number of desired points

%%
% Using |bezierCurve| to get the curve coordinates for n points
V=bezierCurve(p,n);

%%

cFigure; hold on;
hp1=plotV(p,'k.-','lineWidth',lineWidth1,'MarkerSize',markerSize1);
hp2=plotV(V,'r.-','lineWidth',lineWidth2,'MarkerSize',markerSize2);
legend([hp1,hp2],{'Control points','Curve'},'Location','NorthOutside');
axisGeom(gca,fontSize); box on; grid on;
gdrawnow;

%% Example 2: 
c=(4/3)*tan(pi/8);
p=[0 1; c 1; 1 c; 1 0;]; %Control points
n=25; %Number of desired points

%%
% Using |bezierCurve| to get the curve coordinates for n points
V=bezierCurve(p,n);

%%

cFigure; hold on;
hp1=plotV(p,'k.-','lineWidth',lineWidth1,'MarkerSize',markerSize1);
hp2=plotV(V,'r.-','lineWidth',lineWidth2,'MarkerSize',markerSize2);
legend([hp1,hp2],{'Control points','Curve'},'Location','NorthOutside');
view(2); axis tight; axis equal; box on; grid on;
set(gca,'FontSize',fontSize);
gdrawnow;

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
