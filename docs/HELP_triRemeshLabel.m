%% triRemeshLabel
% Below is a demonstration of the features of the |triRemeshLabel| function

%% Syntax
% |[Fb,Vb,Cb]=triRemeshLabel(F,V,pointSpacing);|

%% Description
% Use |triRemeshLabel| to resample and label triangulated surface data.

%% Examples

clear; close all; clc;

%%
% Set path and file name

%Set main folder
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','STL'); 

%File name for STL
fileName=fullfile(pathName,'propeller.stl'); 

%% Import STL surface model
[stlStruct] = import_STL(fileName);

% Access the data from the STL struct
F=stlStruct.solidFaces{1}; %Faces
V=stlStruct.solidVertices{1}; %Vertices

% Merging nodes
[F,V]=mergeVertices(F,V);

%%

pointSpacing=6; 
[Fb,Vb,Cb]=triRemeshLabel(F,V,pointSpacing);

%%
cFigure; 
subplot(1,2,1); hold on; 
title('Original mesh');
gpatch(F,V,'kw','k',1);
axisGeom;
view(-36,-70);
camlight headlight;

subplot(1,2,2); hold on; 
title('Resampled and labelled mesh');
gpatch(Fb,Vb,Cb,'k',1);
axisGeom;
view(-36,-70);
camlight headlight;
colormap(gjet(250)); icolorbar;

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
