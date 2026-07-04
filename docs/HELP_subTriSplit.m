%% subTriSplit
% Below is a demonstration of the features of the |subTriSplit| function

%% Syntax
% |[Fs,Vs]=subTriSplit(F,V,n);|

%% Description
% The |subTriSplit| function enables refinement of triangulated data by
% linearly splitting each triangle into 4 new triangles with each
% iteration. 

%% Examples

clear; close all; clc;

%% 
% Plot Settings
fontSize=15;
faceAlpha=1;
edgeColor=0.2*ones(1,3);
edgeWidth=1.5;
markerSize=35; 
markerSize2=20; 

%% Refining a triangulation through linear splitting

%%
% Example data, an icosahedron
[V,F]=platonic_solid(4);

%%
% Refine using |subTriSplit|

n=3; % Number of split iterations
[Fn,Vn]=subTriSplit(F,V,n); %Refine by linear splitting

%%
cFigure; hold on;
hp1=gpatch(F,V,'none','r',0,4);
% patchNormPlot(F,V,[],[],'r');
hp2=gpatch(Fn,Vn,'w','b',1,1);
% patchNormPlot(Fn,Vn,[],[],'b');
legend([hp1 hp2],{'Original',['Split ',num2str(n),' times' ]});
axisGeom; camlight headlight; 
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
