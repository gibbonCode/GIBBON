%% import_obj
% Below is a demonstration of the features of the |import_obj| function

%%
clear; close all; clc;

%% Syntax
% |objStruct=import_obj(fileName,fullMode);|

%% Description 
% This function imports the OBJ geometry as well as the texture/material
% data when requested. 

%%
clear; close all; clc;

%% Examples
%

%% 
% File name for OBJ file

defaultFolder = fileparts(fileparts(mfilename('fullpath')));
loadPath=fullfile(defaultFolder,'data','OBJ');

testCase=1;
switch testCase
    case 1
        modelName='gibbon.obj';
    case 2
        modelName='test.obj';
    case 3
        modelName='lego_figure.obj';
end
fileName=fullfile(loadPath,modelName);

%% Example 1: Import only the geometry

% Import an OBJ file
objImportOptions.fullMode=0;
objStruct=import_obj(fileName,objImportOptions);
F=objStruct.F;
V=objStruct.V;

%%

cFigure;
gpatch(F,V,'w','k');
axisGeom; camlight headlight;
drawnow;

%% Example 2: Import only the geometry as well as texture data

% Import an OBJ file
objImportOptions.fullMode=1;
objStruct=import_obj(fileName,objImportOptions);
F=objStruct.F;
V=objStruct.V;
C=objStruct.C;
F_uv=objStruct.F_uv;
ij_M=objStruct.ij_M;
m=objStruct.m;

%%

cFigure;
gpatch(F,V,C,'none');
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
