%% triMeshEquilateral
% Below is a demonstration of the features of the |triMeshEquilateral| function

%%
clear; close all; clc;

%%
% PLOT SETTINGS
font_size=20;
cmap=gray(250);
falpha=1;
patch_types={'sx','sy','sz','v'};
ptype=3;
no_slices=4;
mark_siz1=25;
mark_siz2=25;
mark_siz3=15;
line_width1=2;
F_alpha1=1;
F_alpha2=0.8;

%% 
% Control parameters

%Desired mesh point spacing
pointSpacing=1;

%Mesh region extrema
maxV=[5 6];
minV=[-5 -5];

%% CREATING AN EQUILATERAL TRIANGLE MESH

[F,V]=triMeshEquilateral(minV,maxV,pointSpacing);

%%
% Plottting mesh
hf1=cFigure;
title('Equilateral mesh','FontSize',font_size);
xlabel('X','FontSize',font_size);ylabel('Y','FontSize',font_size); zlabel('Z','FontSize',font_size);
hold on;
hpm=patch('Faces',F,'Vertices',V,'EdgeColor','k','FaceColor','g','FaceAlpha',F_alpha1,'lineWidth',line_width1);
% [hp]=patchNormPlot(F,V,1);
colormap autumn; 
axis equal; view(2); axis tight;  grid on; 
set(gca,'FontSize',font_size);
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
