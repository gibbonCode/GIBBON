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
 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
