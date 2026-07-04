%% patchBoundary
% Below is a demonstration of the features of the |patchBoundary| function

%%
clear; close all; clc;

%% Syntax
% |[Eb,E,indBoundary]=patchBoundary(F);|

%% Description 
% This function provides the boundary edges for the input 

%% Examples 
% 

%% Example 1: Get the boundary edges for a mesh

%%
% Creating Example geometry

%Boundary 1
ns=150;
t=linspace(0,2*pi,ns+1);
t=t(1:end-1);
r=6+2.*sin(5*t);
[x,y] = pol2cart(t,r);
V1=[x(:) y(:)];

%Boundary 2
[x,y] = pol2cart(t,ones(size(t)));
V2=[x(:) y(:)+4];

%Boundary 3
[x,y] = pol2cart(t,2*ones(size(t)));
V3=[x(:) y(:)-0.5];

%Defining a region
regionCell={V1,V2,V3}; %A region between V1 and V2 (V2 forms a hole inside V1)
pointSpacing=1; %Desired point spacing
resampleCurveOpt=1; %Option to turn on/off resampling of input boundary curves

[F,V]=regionTriMesh2D(regionCell,pointSpacing,resampleCurveOpt,0);

%%
% Get boundary edges 

Eb=patchBoundary(F);

%%
% Visualize edges

cFigure; 
hp1=gpatch(F,V,'kw','k',1,2);
hp2=gpatch(Eb,V,'none','r',1,3);
legend([hp1 hp2],{'Mesh','Boundary edges'})
axisGeom; view(2);
drawnow;

%% Example 2: Get the boundary edges for a multi mesh type cell
%

[V,F]=patch_dual(V,F);

%%
% Get boundary edges 

Eb=patchBoundary(F);

%%
% Visualize edges

cFigure; 
hp1=gpatch(F,V,'kw','k',1,2);
hp2=gpatch(Eb,V,'none','r',1,3);
legend([hp1(1) hp2],{'Mesh','Boundary edges'})
axisGeom; view(2);
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
