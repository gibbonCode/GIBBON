%% importFEBio_logfile
% Below is a demonstration of the features of the |importFEBio_logfile| function

%%
clear; close all; clc;

%% Syntax
% |dataStruct=importFEBio_logfile(fileNameImport,addZeroInitialOpt,removeIndicesOpt)|

%% Description 
% This function imports FEBio output log files such as for node outputs
% (e.g. displacement) and element outputs (e.g. stress). The inputs include
% the file name to import as well as the following optional inputs: 
% 
% * addZeroInitialOpt => Set to 0 or 1 to import as is, or to add a zero
% initial state. Note this functionality is not obsolute, since FEBio can
% now import an initial state (which may be non-zero e.g. for stretch
% data). 
% * removeIndicesOpt => set to 0 or 1 to keep or remove the index column
% from the data 

%% Examples 
% 

%%

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath'))); 
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress sigma_z

%%
% Check for existing output files in the temp folder. Run febio demo 1 if not found. 

if ~exist(febioLogFileName_disp,'file')
    DEMO_febio_0001_cube_uniaxial
end

%% Importing nodal displacements from a log file
% Importing the output data into a data structure: 

addZeroInitialOpt=0;
removeIndicesOpt=1;
dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),addZeroInitialOpt,removeIndicesOpt)

%%
% The data structure contains the time vector and the data. For k
% dimensional data, n nodes (or elements), and m time steps the data is nxkxm+1 (or m
% when the initial state is not requested or not added by FEBio). 

%Access data
U_disp=dataStruct.data; %Displacement data 
timeVec=dataStruct.time; %Time vector

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
