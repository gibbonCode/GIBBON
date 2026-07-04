%% patchFeatureDetect
% Below is a demonstration of the features of the |patchFeatureDetect| function

%%
clear; close all; clc;

%% Syntax
% |[G]=patchFeatureDetect(F,V,a);|

%% Description 
% This function detects surface features or groups for the input patch data
% defined by the faces F and the vertices V. Patch elements for which the
% dihedral angle is lower than a are grouped together. 
% This function is useful for detecting sets of faces from for instance
% importance CAD geometry, e.g. top faces, side faces etc can be
% automatically detected provided that their boundaries have a significant
% enough angle change. 

%%
% Plot settings
fontSize=15;

%% 
% Create example data 

testCase=4;
switch testCase
    case 1
        boxDim=[4 5 6]; %Width in each direction
        pointSpacing=1; %Desired point spacing
        [F,V]=triBox(boxDim,pointSpacing);
    case 2
        boxDim=[4 5 6]; %Width in each direction
        boxEl=[4 5 6]; %Desired number of elements
        [F,V]=quadBox(boxDim,boxEl);
    case 3
        defaultFolder = fileparts(fileparts(mfilename('fullpath')));
        pathName=fullfile(defaultFolder,'data','libSurf');
        dataStruct=load(fullfile(pathName,'enginePart_p1.mat'));        
        F=dataStruct.F;
        V=dataStruct.V;
    case 4
        defaultFolder = fileparts(fileparts(mfilename('fullpath')));
        pathName=fullfile(defaultFolder,'data','libSurf');
        dataStruct=load(fullfile(pathName,'sprocket.mat'));
        F=dataStruct.F;
        V=dataStruct.V;        
end

%%
% Use |patchFeatureDetect| to detect features in patch

%Angular threshold in radians
a=(45/180)*pi; 

%Detect surface features
G=patchFeatureDetect(F,V,a); 

%%

% Plotting model
cFigure; 
subplot(1,2,1); hold on;
title('Input surface','FontSize',fontSize);
gpatch(F,V,'w','k',1);
axisGeom(gca,fontSize); camlight headlight; 

subplot(1,2,2); hold on;
title('Surface feature labels','FontSize',fontSize);
gpatch(F,V,G,'k',1);
axisGeom(gca,fontSize); camlight headlight; 
colormap gjet; icolorbar; 
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
