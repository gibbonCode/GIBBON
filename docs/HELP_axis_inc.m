%% axis_inc
% Below is a demonstration of the features of the |axis_inc| function

%%
clear; close all; clc;

%% Syntax
% |h=axis_inc(s);|

%% Description 
% This function widens the axis limits for the axis with the handle |hs|
% using the scale factor |s| (e.g. if s=2 the axis width is doubled). The
% scaling parameter |s| can be a single scalar or a vector such that each
% axis direction is scaled differently.  
% 
%% Examples 
% 

%%
% Plot settings
fontSize=25;
lineWidth=3; 

%% Scaling limits on a 2D plot

scaleFactor=2; 

cFigure; 
h1=subplot(1,2,1); hold on; 
title('Original axis limits');
plot([-1 1 1 -1],[-1 -1 1 1],'b-','LineWidth',lineWidth);
axis tight; axis equal; view(2); box on; grid on;
set(h1,'FontSize',fontSize);

h2=subplot(1,2,2); hold on; 
title(['Axis limits expanded using scale factor of ',sprintf('%.2f',scaleFactor)]);
plot([-1 1 1 -1],[-1 -1 1 1],'b-','LineWidth',lineWidth);
axis tight; axis equal; view(2); box on; grid on;
axis_inc(scaleFactor,h2);
set(h2,'FontSize',fontSize);
drawnow; 

%% Scaling limits on a 3D plot

[F,V]=graphicsModels(8);

scaleFactor=[2 4 2]; 

cFigure; 
h1=subplot(1,2,1); hold on; 
title('Original axis limits');
gpatch(F,V,'bw','none');
axisGeom(h1,fontSize);
set(h1,'FontSize',fontSize);
camlight headlight;

h2=subplot(1,2,2); hold on; 
title(['Axis limits expanded using scale factor of ',sprintf('%.2f, %.2f, %.2f',scaleFactor)]);
gpatch(F,V,'bw','none');
axisGeom(h1,fontSize);
axis_inc(scaleFactor,h2);
set(h2,'FontSize',fontSize);
camlight headlight;
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
