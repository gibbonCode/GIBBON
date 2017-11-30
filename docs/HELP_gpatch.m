%% gpatch
% Below is a demonstration of the features of the |gpatch| function

%%
clear; close all; clc;

%% Syntax
% |[hp]=gpatch(F,V,C,CE,A,L);

%% Description
% This function is a short-hand version of the |patch| command. The inputs
% for |gpatch| are the faces (F), the vertices (V), the color description
% (C), the edge color description CE, the transparancy (A), and the edge
% width (L). 
% The color description C can be: 
% 1) A string such as 'g' for green
% 2) A triplet of RGD values e.g. [1 0 0] is blue
% 3) A nx1 or a mx1 array of colormapped colors (where n=size(F,1) or m=size(V,1))
% 4) (simiarl to 3) A nx3 or a mx3 RGB color value array for the faces or vertices respectively. 

%% Examples

%%
% Create example mesh data

[F,V,~]=geoSphere(2,1); %Faces and vertices
CV=V(:,1); %Color information for vertices
CF=vertexToFaceMeasure(F,CV); %Color information for faces
CF_rgb=abs(vertexToFaceMeasure(F,V)); %Color information for faces

%% Example: Introduction to using |gpatch| for mesh visualization
% The below visualization show the syntax require using |patch| and
% |gpatch|. Essentially |gpatch| is just a short-hand version of |patch|
% allowing for quick and easy visualization using patch graphics. 

%%
% Using patch graphics in MATLAB see documentation on |patch| for more information

cFigure; 
suptitle('Using patch')
subplot(2,2,1); 
patch('Faces',F,'Vertices',V,'FaceColor','r','EdgeColor','g','FaceAlpha',0.5);
axisGeom; 

subplot(2,2,2);
patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',CF,'EdgeColor','k','FaceAlpha',1);
colormap gjet;
axisGeom; 

subplot(2,2,3);
patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',CV);
colormap gjet; shading interp;
axisGeom; 

subplot(2,2,4);
patch('Faces',F,'Vertices',V,'FaceColor','flat','FaceVertexCData',CF_rgb,'EdgeColor','k','FaceAlpha',1,'LineWidth',3);
axisGeom; 

drawnow; 

%%
% Using |gpatch| shorthand alternative to |patch|

cFigure; 
suptitle('Using gpatch')
subplot(2,2,1); 
gpatch(F,V,'r','g',0.5);
axisGeom; 

subplot(2,2,2);
gpatch(F,V,CF);
colormap gjet;
axisGeom; 

subplot(2,2,3);
gpatch(F,V,CV);
colormap gjet; shading interp;
axisGeom; 

subplot(2,2,4);
gpatch(F,V,CF_rgb,'k',1,3);
axisGeom; 

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
