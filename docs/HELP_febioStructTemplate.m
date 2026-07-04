%% febioStructTemplate
% Below is a demonstration of the features of the |febioStructTemplate| function

%%
clear; close all; clc;

%% Syntax
% |[domNode]=febioStructTemplate(febio_spec,fileName,optionStruct);|

%% Description
% This function creates a template FEBio input file structure. 

%%
% See also: |febioStruct2xml|

%% Creating a template FEBio input structure
% The created template contains the following sections: module, control,
% globals, loadData, and output
[febio_spec]=febioStructTemplate;

%% Viewing the FEBio structure
% The |febView| command can be used to render the XML form of the febio structure in a figure
% window. 

%%
% NOTE: The figure below does not render in documentation due to a
% MATLAB but (or limitation). The code |[hFig]=febView(domNode);| is
% therefore suppressed.
%%
% |[hFig]=febView(febio_spec);|

%% Export to an FEBio input file
% You can use |febioStruct2xml| to write the xml data to a file AND/OR a
% domNode object (optional output). Leave the fileName variable empty to
% supress file export. 

%Create file name for XML file
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');
fileName=fullfile(savePath,'tempModel.feb');

febioStruct2xml(febio_spec,fileName); %Exporting to file and domNode

%% Viewing a FEBio file
% One can also use |febView| for file febio fil viewing: 

%%
% |febView(fileName);|

%%
% Clear variable for next example
clear febio_spec

%% Overriding defaults

%%
% Create your own settings

% Some custom entries to override default
febio_spec.Control.time_steps=33;
febio_spec.Control.step_size=1/febio_spec.Control.time_steps;

% Some custom fields/entries not part of the default
febio_spec.Control.someField='Some custom thing';
febio_spec.Control.someStruct.someData=pi;

%%
% What the custom set looks like:
febio_spec.Control

%%
% Use |febioStructTemplate| to "fill in the gaps"
[febio_spec]=febioStructTemplate(febio_spec);

%%
% What the set looks like after filling in gaps based on the templete:
febio_spec.Control

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
