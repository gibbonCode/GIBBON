%% polythick
% Below is a demonstration of the features of the |polythick| function

%%
clear; close all; clc;

%% Syntax
% |[Xu Yu Xl Yl]=polythick(x,y,t,nSteps,resampleFactor);|

%% Description 
% This function thickens a polygon defined by x, y. It generates the upper
% (Xu, Yu) and lower (Xl, Yl) coordinates depending on the thickness t. 
%
% The slope is calculated at each coordinate and points are copied upwards
% and downwards (thickening) orthogonal to the local slope.

%% Examples 
% 

%%
% Plot settings
markerSize=15;
lineWidth=2; 
fontSize=25; 

%% Example 1: Thickening a curve

x=linspace(0,2*pi,25)';
y=sin(x);
t=0.25;
nSteps=5;
resampleFactor=1; 

[xu,yu,xl,yl]=polythick(x,y,t,nSteps,resampleFactor);

%%

cFigure;  hold on; 
plot(x,y,'g.-','MarkerSize',markerSize,'LineWidth',lineWidth);
plot(xu,yu,'r.-','MarkerSize',markerSize,'LineWidth',lineWidth);
plot(xl,yl,'b.-','MarkerSize',markerSize,'LineWidth',lineWidth);
axis equal; axis tight; grid on; box on; 
set(gca,'FontSize',fontSize)
drawnow; 

%% Example 2: Thickening a curve with curvature

x=linspace(0,2*pi,25)';
y=5*sin(x);
t=0.25;
nSteps=5;
resampleFactor=5; %Upsample to sample high curvature region better
interpMethod='pchip';

[xu,yu,xl,yl]=polythick(x,y,t,nSteps,resampleFactor,interpMethod);

%%

cFigure;  hold on; 
plot(x,y,'g.-','MarkerSize',markerSize,'LineWidth',lineWidth);
plot(xu,yu,'r.-','MarkerSize',markerSize,'LineWidth',lineWidth);
plot(xl,yl,'b.-','MarkerSize',markerSize,'LineWidth',lineWidth);
axis equal; axis tight; grid on; box on; 
set(gca,'FontSize',fontSize)
drawnow; 

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
