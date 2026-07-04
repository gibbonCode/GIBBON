%% patchAnnotate
% Below is a demonstration of the features of the |patchAnnotate| function

%%
clear; close all; clc;

%% Syntax
% |[hv,hf]=patchAnnotate(F,V,d,varargin);|

%% Description 
% This function annotates patch data in terms of node numbers and face
% numbers. The inputs are the faces array (F), the vertex array (V), and
% the surface offset (d). If d=0 the text annotations are positioned at the
% nodes for the vertices and at the face centres for the faces. If d~=0 the
% the text is places d away from the surface in the direction of the
% surface normal (negative d values are permitted). If d is empty the
% offset defaults to 1/4 the average edge length. 
% Other optional additional inputs are those associated with MATLAB's text
% function (e.g. 'FontSize', 'Color', etc). 
%
% Faces are annotated using a bold font while nodes are annotated using a
% non-bold and italic font. 

%% Examples 

%%
% Plot settings
fontSize=15;
markerSize=20;

%%
% Create example patch data
[F,V,~]=geoSphere(1,1);

%%
% Visualize the patch data and annotate useing |patchAnnotate|

cFigure; hold on;
gpatch(F,V); %Add patch 

[hv,hf]=patchAnnotate(F,V,[],'FontSize',fontSize,'Color','k'); %Annotate

axisGeom;
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
