%% export_mesh
% Below is a demonstration of the features of the |export_mesh| function

%%
clear; close all; clc;

%% Syntax
% |export_mesh(fileName,meshStruct);|

%% Description
% This function exports a .MESH file based on the meshStruct which should
% contain the vertices and elements. 
%
% Elements can be specified as either: 
% meshStruct.elements or meshStruct.E
%
% Vertices can be specified as either: 
% meshStruct.vertices or meshStruct.V or meshStruct.nodes
%
% This function currently only supporte tetrahedral and hexahedral element
% exporting. 
%
% The element type is either inferred from the size of the element matrix
% or can be set using: 
% meshStruct.elementType
% Use meshStruct.elementType='Hexahedral', for 8-noded hexahedral elements
% Use meshStruct.elementType='Tetrahedral', for 4-noded tetrahedral elements
%
% This implementation was created using examples of .MESH files contained
% in: https://github.com/cnr-isti-vclab/HexaLab

%% Examples

%% Exporting a .MESH file
% Creating a solid hexahedral mesh sphere

%%
% Get example data 
testCase=1; 
switch testCase
    case 1 %Cube
        boxDim=[1 1 1];
        boxEl=[2 2 2];
        meshStruct=hexMeshBox(boxDim,boxEl);
    case 2 %Sphere
        %Control settings
        optionStruct.sphereRadius=10;
        optionStruct.coreRadius=5;
        optionStruct.numElementsMantel=5;
        optionStruct.numElementsCore=8;
        optionStruct.makeHollow=0;
        optionStruct.outputStructType=2;
        
        %Creating sphere
        [meshStruct]=hexMeshSphere(optionStruct);
end

%%
% Define file name

%Get toolbox folder paths
defaultFolder = fileparts(fileparts(mfilename('fullpath'))); 
savePath=fullfile(defaultFolder,'data','temp');

%Define file name
fileName=fullfile(savePath,'temp.mesh'); %File name for mesh file

%%
% Export mesh to .MESH file
export_mesh(fileName,meshStruct);

%%

type(fileName)

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
