%% quadBox
% Below is a demonstration of the features of the |quadBox| function

%%
clear; close all; clc;

%%
% PLOT SETTINGS
fig_color='w'; fig_colordef='white';
fontSize=20;
edgeWidth=2;
edgeColor=0.7*ones(1,3);
faceAlpha1=0.5;

%% Creating a quadrilateral mesh of a box

%% 
% Specifying dimensions and number of elements for each direction
boxDim=[4 5 6]; %Width in each direction
boxEl=[3 4 5]; %Number of elements per direction 

%%
% Using |quadBox| to build the patch model

[F,V,faceBoundaryMarker]=quadBox(boxDim,boxEl);

%%
% Plotting model
hf1=figuremax(fig_color,fig_colordef);
title('Box quadrilateral faces and normals','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',faceBoundaryMarker,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
[hp]=patchNormPlot(F,V,mean(boxDim./boxEl));

colormap(jet(6)); colorbar; 
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

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
% Copyright (C) 2017  Kevin Mattheus Moerman
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
