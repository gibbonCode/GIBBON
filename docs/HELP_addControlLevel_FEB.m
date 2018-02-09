%% addControlLevel_FEB
% Below is a demonstration of the features of the |addControlLevel_FEB| function

%%
clear; close all; clc;

%% Syntax
% |[domNode]=addControlLevel_FEB(domNode,FEB_struct);|

%% Description
% This function adds the control information to the input XML
% data (domNode) based on the input febio structure (FEB_struct).

%% Examples

%% Example: Defining control section
% 

%Example data 
numSteps=25;

%Control section
FEB_struct.Control.AnalysisType='static';
FEB_struct.Control.Properties={'time_steps','step_size',...
                               'max_refs','max_ups',...
                               'dtol','etol','rtol','lstol'};

FEB_struct.Control.Values={numSteps,1/numSteps,25,0,0.001,0.01,0,0.9};
FEB_struct.Control.TimeStepperProperties={'dtmin','dtmax','max_retries','opt_iter'};
FEB_struct.Control.TimeStepperValues={(1/(100*numSteps)),1/numSteps,5,10};


%Initialize docNode object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec');

% %Add geometry information (Surfaces)
% domNode=addGeometryLevel_FEB(domNode,FEB_struct);

%Add boundary condition information
domNode=addControlLevel_FEB(domNode,FEB_struct);

%%
%  View example XML string
XML_str = xmlwrite(domNode);
disp(XML_str);

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
% Copyright (C) 2018  Kevin Mattheus Moerman
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
