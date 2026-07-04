%% tri2quad
% Below is a demonstration of the features of the |tri2quad| function

%%
clear; close all; clc;

%% Syntax
% |[Fq,Vq]=tri2quad(Ft,Vt);|

%% Description 
% Tris func
%% Examples 
% 

%%
% Plot settings
fontSize=15;

%%
% Example triangulated surface
[F,V]=stanford_bunny;

%% Example 1: Converting a triangulated surface to a quandrangulated surface

%%
% Convert triangular faces to quadrilateral faces
[Fq,Vq]=tri2quad(F,V);

%%
% Visualisation

cFigure; 
subplot(1,2,1);
title('Triangulation','FontSize',fontSize);
gpatch(F,V,'rw','k');
set(gca,'FontSize',fontSize);
axisGeom;
camlight('headlight'); 

subplot(1,2,2);
title('Quadrangulation','FontSize',fontSize);
gpatch(Fq,Vq,'gw','k');
set(gca,'FontSize',fontSize);
axisGeom;
camlight('headlight'); 

drawnow; 

%% Example 2: Using incentres rather than triangle centres

%%
% Convert triangular faces to quadrilateral faces
centerType=2;
[Fq,Vq]=tri2quad(F,V,centerType);

%%
% Visualisation

cFigure; 
subplot(1,2,1);
title('Triangulation','FontSize',fontSize);
gpatch(F,V,'rw','k');
set(gca,'FontSize',fontSize);
axisGeom;
camlight('headlight'); 

subplot(1,2,2);
title('Quadrangulation','FontSize',fontSize);
gpatch(Fq,Vq,'gw','k');
set(gca,'FontSize',fontSize);
axisGeom;
camlight('headlight'); 

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
