%% groupVertices
% Below is a demonstration of the features of the |groupVertices| function

%%
clear; close all; clc;

%% Syntax
% |[groupIndexVertices,groupIndexFaces]=groupVertices(F,V,waitBarOption);|

%% Description
% This function takes in a description of a tesselation, consisting of an
% array F (faces or elements), and V (vertices). The output is represents
% the group index for each point in groupIndexVertices and the group index
% of each face in groupIndexFaces. Vertices which are connected to the same
% mesh obtain the same group indes. Face indices are derived from the
% vertex indices. 

%% Examples

%% Example 1: Seperating patch data into groups

%%
% Creating test data consisting of seperates sets faces and vertices

numGroups=3; 
F_cell=cell(1,numGroups); 
V_cell=cell(1,numGroups); 
for q=1:1:numGroups
    switch q
        case 1
            [F,V]=geoSphere(3,1);
        case 2
            [F,V]=stanford_bunny('g');
        case 3
            [F,V]=graphicsModels(3);
    end
%     [F,V]=subtri(F,V,1);
    V=V-mean(V,1);
    V=V./max(V(:));
    V(:,1)=V(:,1)+q*2;
    F_cell{q}=F;
    V_cell{q}=V;
end
[F,V]=joinElementSets(F_cell,V_cell);

%%
% Using |groupVertices| to split the vertices into groups
[groupIndexVertices,groupIndexFaces]=groupVertices(F,V,1);

%%
% Visualizing the resulting grouping

cFigure; 
subplot(2,1,1); hold on;
title('Ungrouped')
gpatch(F,V,'kw','none');
axisGeom;
camlight headlight; 

subplot(2,1,2); hold on;
title('Grouped')
gpatch(F,V,'kw','none');
scatterV(V,15,groupIndexVertices,'filled');
axisGeom;
camlight headlight; 
colormap gjet; icolorbar;
drawnow;

%%
% If desired a second output can be requisted which represents the face
% groupings. 

cFigure; 
title('Grouped faces')
gpatch(F,V,groupIndexFaces,'none');
axisGeom;
camlight headlight; 
colormap gjet; icolorbar;
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
