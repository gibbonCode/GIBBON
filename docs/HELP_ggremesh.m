
%% ggremesh
% Below is a demonstration of the features of the |ggremesh| function

%% Syntax
% |[Fn,Vn]=ggremesh(F,V,optionStruct)|

%% Description 
% This function uses the external library Geogram to remesh the input
% triangulation defined by the faces F and the vertices V. In particular
% the code "vorpalite" is used. An additional option structure may be
% provided where users can set particular parameters for Geogram. 
%
% Below the options and defaults are provided: 
% optionStruct.nb_pts=size(V,1); %number of points
% optionStruct.anisotropy=0; %Use anisotropy (~=0) to capture geometry or favour isotropic triangles (=0)
% optionStruct.pre.max_hole_area=100; %Max hole area for pre-processing step
% optionStruct.pre.max_hole_edges=0; %Max number of hole edges for pre-processing step
% optionStruct.post.max_hole_area=100; %Max hole area for post-processing step
% optionStruct.post.max_hole_edges=0; %Max number of hole edges for post-processing step
% optionStruct.disp_on=1; %Turn on/off displaying of Geogram text
%
% Instead of nb_pts users can also specify a pointSpacing to be used
% instead of nb_pts. This is not a Geogram feature but a GIBBON option
% which is translated to the number of points for Geogram remeshing. This
% is and example for a desired point spacing of 4:  
% optionStruct.pointSpacing=4
%
% Geogram website:
% http://alice.loria.fr/index.php/software/4-library/75-geogram.html 
% 
% Geogram license: 
% http://alice.loria.fr/software/geogram/doc/html/geogram_license.html
%
% Lévy B., Bonneel N. (2013) Variational Anisotropic Surface Meshing with
% Voronoi Parallel Linear Enumeration. In: Jiao X., Weill JC. (eds)
% Proceedings of the 21st International Meshing Roundtable. Springer,
% Berlin, Heidelberg. https://doi.org/10.1007/978-3-642-33573-0_21 
% 
% See also: 
% http://alice.loria.fr/publications/papers/2012/Vorpaline_IMR/vorpaline.pdf
% https://www.ljll.math.upmc.fr/hecht/ftp/ff++days/2013/BrunoLevy.pdf

%% Examples 

%%
clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceColor='b';
faceAlpha=1;
edgeColor='k';
edgeWidth=0.5;

%% Example 1: Remeshing a triangulated surface isotropically

%% 
% Get example geometry
[F,V]=graphicsModels(5); % Get surface

%%
% Remesh using ggremesh
[Fn,Vn]=ggremesh(F,V);

%%
% Visualiza patch data

cFigure; 
subplot(1,2,1); hold on;
title('Input mesh');
gpatch(F,V,'w','k');
axisGeom;
view(-75,-36);
camlight headlight; axis off;

subplot(1,2,2); hold on;
title('Geogram remeshed');
gpatch(Fn,Vn,'gw','k',1,1);
axisGeom;
view(-75,-36);
camlight headlight; axis off;

gdrawnow;

%% Example 2: Remeshing a triangulated surface with desired number of points

%% 
% Get example geometry
[F,V]=stanford_bunny; % Get surface

%%
% Remesh using ggremesh

optionStruct1.nb_pts=1500; %Set desired number of points
optionStruct1.disp_on=1;
[Fn,Vn]=ggremesh(F,V,optionStruct1);

%%
% Visualiza patch data

cFigure; 
subplot(1,2,1); hold on;
title('Input mesh');
gpatch(F,V,'w','k');
axisGeom;
camlight headlight; axis off;

subplot(1,2,2); hold on;
title('Geogram remeshed');
gpatch(Fn,Vn,'gw','k',1,1);
axisGeom;
camlight headlight; axis off;

gdrawnow;

%% Example 3: Remeshing a triangulated surface with desired point spacing

%% 
% Get example geometry
[F,V]=graphicsModels(11); % Get surface

%%
% Remesh using ggremesh

optionStruct2.pointSpacing=3; %Set desired point spacing
optionStruct2.disp_on=1; % Turn off command window text display
% optionStruct2.anisotropy=5;
[Fn,Vn]=ggremesh(F,V,optionStruct2);

%%
% Visualiza patch data

cFigure; 
subplot(1,2,1); hold on;
title('Input mesh');
gpatch(F,V,'w','k');
axisGeom;
camlight headlight; 

subplot(1,2,2); hold on;
title('Geogram remeshed');
gpatch(Fn,Vn,'gw','k',1,1);
axisGeom;
camlight headlight; 

gdrawnow;

%% Example 4: Setting pre- and prost-processing settings e.g. to close holes

%% 
% Get example geometry

inputStruct.cylRadius=1;
inputStruct.numRadial=15;
inputStruct.cylHeight=3;
inputStruct.numHeight=11;
inputStruct.meshType='tri';

% Derive patch data for a cylinder
[F,V]=patchcylinder(inputStruct); 

%% 
% Remesh using ggremesh
optionStruct3.nb_pts=size(V,1); %Set desired number of points
optionStruct3.disp_on=0; % Turn off command window text display
optionStruct3.pre.max_hole_area=10; %Max hole area for pre-processing step
optionStruct3.pre.max_hole_edges=20; %Max number of hole edges for pre-processing step
% optionStruct3.post.max_hole_area=10; %Max hole area for post-processing step
% optionStruct3.post.max_hole_edges=20; %Max number of hole edges for post-processing step

[Fn,Vn]=ggremesh(F,V,optionStruct3);

% Visualiza patch data
Eb=patchBoundary(F);
cFigure; 
subplot(1,2,1); hold on;
title('Input mesh with holes');
gpatch(F,V,'w','k');
gpatch(Eb,V,'none','b',1,2);
axisGeom;
camlight headlight; 

subplot(1,2,2); hold on;
title('Geogram remeshed and closed');
gpatch(Fn,Vn,'gw','k',1,1);
axisGeom;
camlight headlight; 

gdrawnow;

%% Example 5: Setting pre- and prost-processing settings e.g. to avoid closure of holes

%% 
% Get example geometry

inputStruct.cylRadius=1;
inputStruct.numRadial=15;
inputStruct.cylHeight=3;
inputStruct.numHeight=11;
inputStruct.meshType='tri';

% Derive patch data for a cylinder
[F,V]=patchcylinder(inputStruct); 

%% 
% Remesh using ggremesh
optionStruct3.nb_pts=size(V,1); %Set desired number of points
optionStruct3.disp_on=1; % Turn off command window text display
optionStruct3.pre.max_hole_area=100; %Max hole area for pre-processing step
optionStruct3.pre.max_hole_edges=0; %Max number of hole edges for pre-processing step

[Fn,Vn]=ggremesh(F,V,optionStruct3);

% Visualiza patch data
Eb=patchBoundary(F);
cFigure; 
subplot(1,2,1); hold on;
title('Input mesh with holes');
gpatch(F,V,'w','k');
gpatch(Eb,V,'none','b',1,2);
axisGeom;
camlight headlight; 

subplot(1,2,2); hold on;
title('Geogram remeshed with holes');
gpatch(Fn,Vn,'gw','k',1,1);
axisGeom;
camlight headlight; 

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
