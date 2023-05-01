%% tri2quad
% Below is a demonstration of the features of the |tri2quad| function

%% Syntax
% |[Fq,Vq]=tri2quad(Ft,Vt);|

%% Description 
% 
%% Examples 
% 
%%
clear; close all; clc;

% Plot settings
fig_color='w'; fig_colordef='white';
fontSize=15;
faceColor='b';
faceAlpha=1;
edgeColor='k';
edgeWidth=1;

%% Example: Converting a triangulated surface to a quandrangulated surface

%%
% Example triangulated surface
[F,V]=stanford_bunny;

%%
% Convert triangular faces to quadrilateral faces
[Fq,Vq]=tri2quad(F,V,1);

%%
% Visualisation

hf=cFigure; 
subplot(1,2,1);
title('Triangulation','FontSize',fontSize);
gpatch(F,V,'rw','r');
set(gca,'FontSize',fontSize);
axisGeom;
camlight('headlight'); 

subplot(1,2,2);
title('Quadrangulation','FontSize',fontSize);
gpatch(Fq,Vq,'gw','g');
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
