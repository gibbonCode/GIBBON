%% patchEdges
% Below is a demonstration of the features of the |patchEdges| function

%%
clear; close all; clc;

%% Syntax
% | [E] = patchEdges(F,uniOpt);|

%% Description 
% Get the nx2 patch array for the input faces (F) and vertices (V). If
% uniOpt=1 then only unique edges are returned e.g. the edges spanning from
% node 1 to node 2 is the same as from node 2 to node 1. If uniOpt=0 then
% such non-unique edges are also inluded. 

%% Examples 
% 

%%
% Example geometry
[F,V]=geoSphere(1,1);

uniOpt=0;
E=patchEdges(F,uniOpt);

cFigure; 
subplot(1,2,1); 
title('Patch geometry')
gpatch(F,V,'rw','r'); 
axisGeom; camlight headlight; 

subplot(1,2,2); 
title('Edges')
gedge(E,V,'r',2);
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
