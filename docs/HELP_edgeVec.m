%% edgeVec
% Below is a demonstration of the features of the |edgeVec| function

%%
clear; close all; clc;

%% Syntax
% |[N,Vp,Nv]=edgeVec(E,V);|

%% Description 
% Computes the edge vectors for the input edges defined by the edge array E
% and the vertices V. The ouput is the allong edge unit vectors N. Other
% optional outputs include the edge vector origins Vp, and the vertex edge
% unit vectors. The latter are an average of both edges connected to the
% vertex. 

%% Examples 
% 

%%
% Plot settings
fontSize=15;
lineWidth1=2; 
lineWidth2=3; 

%%
% Example edge data

%Example mesh
plateDim = [1 1]; 
plateEl  = [2 3];
[F,V]=quadPlate(plateDim,plateEl);

%Get boundary edges as example
E=patchBoundary(F);

%%
% Get edge vectors

[N,Vp,Nv]=edgeVec(E,V);

%%
% Visualization

cFigure; 
subplot(1,2,1); hold on; 
hp1=gpatch(F,V,'kw','w',0.5,lineWidth1);
hp2=gedge(E,V,'k',lineWidth2);
hp3=quiverVec(Vp,N,0.2,'o');
legend([hp1 hp2 hp3],{'Mesh','Input edges','Output edge vectors'},'Location','northoutside');
axisGeom(gca,fontSize); view(2);

subplot(1,2,2); hold on; 
hp1=gpatch(F,V,'kw','w',0.5,lineWidth1);
hp2=gedge(E,V,'k',lineWidth2);
hp3=quiverVec(V,Nv,0.2,'o');
legend([hp1 hp2 hp3],{'Mesh','Input edges','Output edge vectors at vertices'},'Location','northoutside');
axisGeom(gca,fontSize); view(2);

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
