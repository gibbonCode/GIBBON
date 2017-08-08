%% addConstraintsLevel_FEB
% Below is a demonstration of the features of the |addConstraintsLevel_FEB| function

%%
clear; close all; clc;

%% Syntax
% |[domNode]=addConstraintsLevel_FEB(domNode,FEB_struct);|

%% Description
% This function adds the contraints information to the input XML
% data (domNode) based on the input febio structure (FEB_struct).

%% Examples

%% Example: Defining fix/prescribed constraint for a rigid body
% 

%Example data 
bcPrescribeMagnitudes=[pi 1 -5]; %example prescribed displacement magnitudes
rigidMaterialId=2; %Example material Id for rigid body material

%Constraint section
FEB_struct.Constraints{1}.RigidId=rigidMaterialId;

FEB_struct.Constraints{1}.Prescribe{1}.bc='x';
FEB_struct.Constraints{1}.Prescribe{1}.Scale=bcPrescribeMagnitudes(1);
FEB_struct.Constraints{1}.Prescribe{1}.lc=1;

FEB_struct.Constraints{1}.Prescribe{2}.bc='y';
FEB_struct.Constraints{1}.Prescribe{2}.Scale=bcPrescribeMagnitudes(2);
FEB_struct.Constraints{1}.Prescribe{2}.lc=1;

FEB_struct.Constraints{1}.Prescribe{3}.bc='z';
FEB_struct.Constraints{1}.Prescribe{3}.Scale=bcPrescribeMagnitudes(3);
FEB_struct.Constraints{1}.Prescribe{3}.lc=1;

FEB_struct.Constraints{1}.Fix{1}.bc='Rx';
FEB_struct.Constraints{1}.Fix{2}.bc='Ry';
FEB_struct.Constraints{1}.Fix{3}.bc='Rz';

%Initialize docNode object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec');

%Add boundary condition information
domNode=addConstraintsLevel_FEB(domNode,FEB_struct);

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
% ********** _license boilerplate_ **********
% 
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
