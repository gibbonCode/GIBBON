%% tri2rhombi
% Below is a demonstration of the features of the |tri2rhombi| function

%%
clear; close all; clc;

%% Syntax
% |[Ft,Vt,Ct]=tri2rhombi(F,V,C);|

%% Description
% This function refines the surface region defined by the logic logicFaces,
% belonging to the surface given by the faces F and vertices V. The region
% is refined by: 1) taking the dual tesselation, 2) retriangulating the
% dual mesh to include the original points.
%
% The input consists of:
% F:          the faces
% V:          the vertices
% logicFaces: A logic for the faces requiring refinement
% C:          color data on either the faces or the vertices
% 
% The ouput can consist of: 
% FT:     Faces
% VT:     Vertices
% C_type: Color/label for triangle type, i.e. original (1), refined (2), or
%         boundary (3) 
% indIni: Indices for original points
% C_new:  New color data for faces or vertices


%% Examples

%%
% Plot settings
fontSize=15;
cMap=gjet(4);
faceAlpha=0.5;
plotColor1=cMap(1,:);
plotColor2=cMap(2,:);
plotColor3=cMap(4,:);
faceColor1=0.5*ones(1,3);
edgeWidth=1;
markerSize=25;

%% Example: Converting triangles to rhombi (quads)

%%
% Building example geometry
[F,V]=graphicsModels(9);

%%    
% Refine surface region using tri2rhombi
[Ft,Vt]=tri2rhombi(F,V);

%%

%Plotting results
cFigure;
subplot(1,2,1); hold on;
title('Original','FontSize',fontSize);
gpatch(F,V,'rw','k',1,edgeWidth);
axisGeom(gca,fontSize);
camlight headlight;
view(0,0); zoom(2);
axis off; 

subplot(1,2,2); hold on;
title('Rhombi','FontSize',fontSize);
gpatch(Ft,Vt,'gw','k',1,edgeWidth);
axisGeom(gca,fontSize);
camlight headlight;
view(0,0); zoom(2);
axis off;
drawnow;

%% Example: Converting triangles to rhombi on a non-closed surface

%%
% Building example geometry

%Defining geodesic dome
r=1; %sphere radius
n=1; %Refinements
[F,V,~]=hemiSphereMesh(n,r,0);

%%    
% Refine surface region using tri2rhombi
[Ft,Vt]=tri2rhombi(F,V);

%%

%Plotting results
cFigure;
subplot(1,2,1); hold on;
title('Original','FontSize',fontSize);
gpatch(F,V,'rw','k',1,edgeWidth);
axisGeom(gca,fontSize);
camlight headlight;
axis off; 

subplot(1,2,2); hold on;
title('Rhombi','FontSize',fontSize);
gpatch(Ft,Vt,'gw','k',1,edgeWidth);
axisGeom(gca,fontSize);
camlight headlight;
axis off;
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
