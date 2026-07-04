%% import_INP
% Below is a demonstration of the features of the |import_INP| function

%%
clear; close all; clc;

%%
% Plot settings
faceAlpha=0.5;
fontSize=25; 

%% Specifying file location

%Set main folder
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','INP'); 

%Set model
testCase=3;
switch testCase
    case 1 %Triangular elements
        fileNameEnd='example_TRI.inp';
        numberNodesElement=3;
    case 2 %Quad elements
        fileNameEnd='example_QUAD.inp';
        numberNodesElement=4;
    case 3 %Tetrahedral elements
        fileNameEnd='example_TET.inp';
        numberNodesElement=4;       
    case 4 %Hexahedral elements
        fileNameEnd='example_HEX.inp';
        numberNodesElement=8;
end
fileName=fullfile(pathName,fileNameEnd); 

%% Importing the INP file
logicRenumberOption=1; 
[elementStruct,nodeStruct]=import_INP(fileName,numberNodesElement,logicRenumberOption);

%%
% Content:
elementStruct
nodeStruct

V=nodeStruct.N; %The nodes
E=elementStruct.E; %The elements

%% 
% Displaying the model
%Get patch data for plotting
if ~isempty(strfind(elementStruct.E_type,'S4R')) || ~isempty(strfind(elementStruct.E_type,'STRI3')); %quad or tri elements
    F=E; %elements already describe faces
else %hex or tet elements
    [F,~]=element2patch(E,[]);    
end

cFigure;
title('INP imported model','fontSize',fontSize);
gpatch(F,V,'bw','b');
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
