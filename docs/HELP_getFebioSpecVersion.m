%% getFebioSpecVersion
% Below is a demonstration of the features of the |getFebioSpecVersion| function

%%
clear; close all; clc;

%% Syntax
% |[febio_spec]=getFebioSpecVersion(febXML);|
% |[febio_spec]=getFebioSpecVersion(febFileName);|

%% Description 
% Returns the febio_spec version used in a particular FEB file specified by
% the input variable febXML. febXML may be of the XML class or specify a
% file name for a feb file. 

%% Examples 
% 

%%
% set up path and file names
defaultFolder = fileparts(fileparts(mfilename('fullpath'))); 
loadPath=fullfile(defaultFolder,'data','FEB');
febioFebFileNamePart='tempModel';

%%
% Get the febio spec from a FEB file

febioFebFileName=fullfile(loadPath,[febioFebFileNamePart,'_4p0.feb']); %FEB file name
getFebioSpecVersion(febioFebFileName)

febioFebFileName=fullfile(loadPath,[febioFebFileNamePart,'_3p0.feb']); %FEB file name
getFebioSpecVersion(febioFebFileName)

febioFebFileName=fullfile(loadPath,[febioFebFileNamePart,'_2p0.feb']); %FEB file name
getFebioSpecVersion(febioFebFileName)

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
