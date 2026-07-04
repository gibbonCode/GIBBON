%% dicom3Dpar
% Below is a demonstration of the features of the |dicom3Dpar| function

%%
clear; close all; clc;

%% Syntax
% |[varargout]=dicom3Dpar(D);|

%% Description 
% This function retrieves the voxel size, position, and orientation
% information from a DICOM file.  

%% Examples 
% 

%% 
% Path name for dicom files

defaultFolder = fileparts(fileparts(mfilename('fullpath'))); %Set main folder
pathName=fullfile(defaultFolder,'data','DICOM','0001_human_calf');
fileName=fullfile(pathName,'00001.dcm');

%%
% Get dicom information
dcmInfo=dicominfo(fileName);

%%
% Get 3D/geometry information, e.g. voxel size, position, and orientation
% information. 

[v,OR,r,c]=dicom3Dpar(dcmInfo)

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
