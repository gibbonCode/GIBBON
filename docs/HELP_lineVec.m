%% lineVec
% Below is a demonstration of the features of the |lineVec| function

%%
clear; close all; clc;

%% Syntax
% |[h]=lineVec(P,V,vecSize,colorSpec,quiverStyleOpt,alphaLevel)|

%% Description
% This function visualizes vector or direction data as lines. Lines can be
% colormapped. 
% The inputs should include the vector origins P, and the vectors V. 
% The following optional inputs can be provided: 
%
% vecSize: The vector size either a single number of [size(V,1) 1]
% colorSpec: string color options e.g. 'k' for black or [size(V,1) 1] or
% [size(V,1) 3] color data 
% quiverStyleOpt: Style option
% lineWidth: Line width
% alphaLevel: Transparency level
%
% The quiver style option can be set to:
% 1: Depart from origin
% 2: Arrive at origin
% 3: Pass through origin
% 4: Two-sided
            
%% Examples

%%
% Plot settings
fontSize=20;

%% Example: basic use for normal direction vector plotting

%%
% Create example data 

[F,V]=geoSphere(2,2);
[N,VN]=patchNormal(F,V);

%%
% Visualizing lines/vectors

cFigure; hold on;
gpatch(F,V,'w','k');
h=lineVec(VN,N);
axisGeom(gca,fontSize); 
camlight headlight; 
gdrawnow; 

%% Example: Specifying full input set

vecSize=0.25; %Vector size
colorSpec=VN(:,3); %Color data 
quiverStyleOpt=3; %Style option
lineWidth=3; %Line width
alphaLevel=0.8; %Transparency level

%%
% Plotting model
cFigure; hold on;
gpatch(F,V,'w','none',0.5);
h=lineVec(VN,N,vecSize,colorSpec,quiverStyleOpt,lineWidth,alphaLevel);
axisGeom(gca,fontSize); 
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
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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
