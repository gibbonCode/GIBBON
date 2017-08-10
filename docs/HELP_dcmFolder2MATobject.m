%% dcmFolder2MATobject
% Below is a demonstration of the features of the |dcmFolder2MATobject| function

%% Syntax
% |dcmFolder2MATobject(PathName,MaxVarSize);|

%% Description 
% The |dcmFolder2MATobject| function converts DICOM data to a MATLAB mat (object and
% or) file. 

%% Examples

clear; close all; clc;

%%
% Plot settings
fig_color='w'; fig_colordef='white';

%% Example: CONVERTING DICOM IMAGE DATA TO A MAT OBJECT
% Below some example code is shown to convert all DICOM files inside a
% folder (including its subfolders) to the IMDAT format. The function
% |dcmFolder2MATobject| converts the DICOM data to a matlab MAT object and
% saves it under the name IMDAT.mat inside a subfolder called IMDAT.
% A |waitbar| appears showing the process of the data conversion for the
% DICOM information and image data. Multiple types of image data (e.g.
% phase, real, imaginary, magnitude data) is stored seperately. Also
% several DICOM info fields are harvested and stored.
%
% The IMDAT.mat object contains the following:
%
% IMDAT_struct =
%
%                 G: [1x1 struct] %The geometry parameters
%         ImageSize: [128 128 17 20] % The image size
%     ImageTypesUni: {'ORIGINAL\PRIMARY\M_FFE\M\FFE'} % The image type or types
%            type_1: [4-D uint16] % The image data matrix
%       type_1_info: [1x340 struct] % The harvested DICOM information
%
% The geometry set G contains:
% G =
%
%      v: [3x1 double] %The voxel size
%     OR: [3x1 double] %The location of the origin
%      r: [3x1 double] %Direction vector for rows
%      c: [3x1 double] %Direction vector for Columns

%%

defaultFolder = fileparts(fileparts(mfilename('fullpath'))); %Set main folder
pathName=fullfile(defaultFolder,'data','DICOM','KNEE_UTE');

%Get all subfolders
if ispc
    pathNames = regexp(genpath(pathName),[filesep,';'], 'split');
elseif isunix
    pathNames = regexp(genpath(pathName),':', 'split');
else
    pathNames = regexp(genpath(pathName),[filesep,';'], 'split');
end

pathNames=pathNames(1:end-1)';
numberOfFolders=numel(pathNames);

%Converting DICOM data to IMDAT format in all subfolders
for q=1:1:numberOfFolders
    pathNameSub=pathNames{q}; %Current path name
    if isempty(strfind(pathNameSub,[filesep,'IMDAT'])); %if the IMDAT directory does not exist yet
        try
            dcmFolder2MATobject(pathNameSub,[]);%Get DICOM data
        catch exception
            warning([exception.message,'-> Analysis skipped for: ',pathName]);
        end
    end
end

%% Example: LOADING OR HANDLING THE MAT OBJECT
% Here is an example for loading in the entire data structure

loadName=fullfile(pathName,'IMDAT','IMDAT.mat');
IMDAT_struct=load(loadName);

%%
% Indexing into the MAT object to avoid loading entire structure
% In somecases it is not desirable to load in the entire data set but only
% say a certain slice. In this case the MAT object allows for indexing as
% shows below. See also the help documentation for |matfile|
% Although this type of indexing can be slow it does allow one to only
% select a subset of the data which in some cases helps to save memory

matObj = matfile(loadName);
G = matObj.G;
M= matObj.type_1;

%% 
% Viewing the image data
sv3(M);

%% 
% Viewing the image data using |ind2patch|
% Alternatively the image data can be viewed using the |ind2patch|
% function. See the associated help for more information.

%%
%
% <<gibbVerySmall.gif>>
%
% _*GIBBON*_
% <www.gibboncode.org>
%
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% <-- GIBBON footer text --> 
% 
%     GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
%     image segmentation, image-based modeling, meshing, and finite element
%     analysis.
% 
%     Copyright (C) 2017  Kevin Mattheus Moerman
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
