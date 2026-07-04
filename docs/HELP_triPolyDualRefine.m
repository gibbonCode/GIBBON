%% triPolyDualRefine
% Below is a demonstration of the features of the |triPolyDualRefine| function

%%
clear; close all; clc;

%% Syntax
% |[Ft,Vt,Ct,indIni]=triPolyDualRefine(F,V);|

%% Description
% This function refines a triangulated polyhedron by first creating the
% dual tesselation and by then triangulating the dual tesselation by
% reintroducing the original point set near the centres of the dual faces. 

%% Examples

%%
% Plot settings

fontSize=15;
faceAlpha=1;
edgeColor=0.*ones(1,3);
edgeWidth=1;
markerSize=15; 

%% Example: Subtriangulating a triangulated surface
% Building example geometry
[F,V]=regionTriMesh2D({[-1 -1; -1 1; 1 1; 1 -1]},0.25,1,0); 
V(:,3)=0;

%%

[Fq,Vq,Cq,indIni]=triPolyDualRefine(F,V);

%%

cMap=gjet(size(Fq,1));
cMap=cMap(randperm(size(cMap,1)),:);

cFigure; 
subplot(1,2,1); hold on; 
title('Original','FontSize',fontSize); 
gpatch(F,V,'gw','k');
axisGeom(gca,fontSize);
view(2);

subplot(1,2,2); hold on; 
title('Refined','FontSize',fontSize); 
gpatch(Fq,Vq,Cq,'k');
plotV(Vq(indIni,:),'k.','MarkerSize',markerSize);
axisGeom(gca,fontSize);
colormap(cMap);
view(2);
drawnow; 

%% Example: Subtriangulating a closed polyhedron (sphere)
% Building example geometry

%Defining geodesic dome
r=1; %sphere radius
n=2; %Refinements   
[F,V,~]=geoSphere(n,r);

%%

[Fq,Vq,Cq,indIni]=triPolyDualRefine(F,V);

%%

cMap=gjet(size(Fq,1));
cMap=cMap(randperm(size(cMap,1)),:);

cFigure; 
subplot(1,2,1); hold on; 
title('Original','FontSize',fontSize); 
gpatch(F,V,'gw','k');
axisGeom(gca,fontSize);
camlight headlight;

subplot(1,2,2); hold on; 
title('Refined','FontSize',fontSize); 
gpatch(Fq,Vq,Cq,'k');
plotV(Vq(indIni,:),'k.','MarkerSize',markerSize);
axisGeom(gca,fontSize);
camlight headlight;
colormap(cMap);
drawnow; 

%% Example: Subtriangulating a closed polyhedron (dinosaur)
% Building example geometry

[F,V]=parasaurolophus;

%%

[Fq,Vq,Cq,indIni]=triPolyDualRefine(F,V);

%%

cMap=gjet(size(Fq,1));
cMap=cMap(randperm(size(cMap,1)),:);

cFigure; 
subplot(1,2,1); hold on; 
title('Original','FontSize',fontSize); 
gpatch(F,V,'gw','k');
axisGeom(gca,fontSize);
camlight headlight;

subplot(1,2,2); hold on; 
title('Refined','FontSize',fontSize); 
gpatch(Fq,Vq,Cq,'k');
plotV(Vq(indIni,:),'k.','MarkerSize',markerSize);
axisGeom(gca,fontSize);
camlight headlight;
colormap(cMap);
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
