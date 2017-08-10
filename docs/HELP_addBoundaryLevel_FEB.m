%% addBoundaryLevel_FEB
% Below is a demonstration of the features of the |addBoundaryLevel_FEB| function

%%
clear; close all; clc;

%% Syntax
% |[domNode]=addBoundaryLevel_FEB(domNode,FEB_struct);|

%% Description
% This function adds the boundary condition information to the input XML
% data (domNode) based on the input febio structure (FEB_struct).

%% Examples

%% Example: Prescribed dislacement using node list
% This method of defining a boundary conditions allows one to vary the
% nodeScale of each member of the set list 

%Prescribed displacement in x direction
FEB_struct.Boundary.Prescribe{1}.Set=[1:10]';
FEB_struct.Boundary.Prescribe{1}.bc='x';
FEB_struct.Boundary.Prescribe{1}.lc=1;
FEB_struct.Boundary.Prescribe{1}.nodeScale=rand(10,1);

%Initialize docNode object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec');

%Add boundary condition information
domNode=addBoundaryLevel_FEB(domNode,FEB_struct);

%%
%  View example XML string
XML_str = xmlwrite(domNode);
disp(XML_str);

%% Example: Prescribed dislacement using NodeSet
% This method of defining a boundary conditions will set a same scale for
% all members of NodeSet

clear FEB_struct; %Clear variable from last example

%Defining node set
FEB_struct.Geometry.NodeSet{1}.Set=11:20;
FEB_struct.Geometry.NodeSet{1}.Name='bcPrescribeList';

%Prescribed displacement in x direction
FEB_struct.Boundary.Prescribe{1}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
FEB_struct.Boundary.Prescribe{1}.bc='x';
FEB_struct.Boundary.Prescribe{1}.lc=1;
FEB_struct.Boundary.Prescribe{1}.Scale=0.5;

%Initialize docNode object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec');

% %Add geometry information (NodeSet)
% domNode=addGeometryLevel_FEB(domNode,FEB_struct);

%Add boundary condition information
domNode=addBoundaryLevel_FEB(domNode,FEB_struct);

%%
%  View example XML string
XML_str = xmlwrite(domNode);
disp(XML_str);

%% Example: Fixed conditions using list

clear FEB_struct; %Clear variable from last example

FEB_struct.Boundary.Fix{1}.bc='x';
FEB_struct.Boundary.Fix{1}.Set=1:10;

%Initialize docNode object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec');

%Add boundary condition information
domNode=addBoundaryLevel_FEB(domNode,FEB_struct);

%%
% View example XML string
XML_str = xmlwrite(domNode);
disp(XML_str);

%% Example: Fixed conditions using NodeSet

clear FEB_struct; %Clear variable from last example

%Defining node set
FEB_struct.Geometry.NodeSet{1}.Set=11:20;
FEB_struct.Geometry.NodeSet{1}.Name='bcFixList';

FEB_struct.Boundary.Fix{1}.bc='x';
FEB_struct.Boundary.Fix{1}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;

%Initialize docNode object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec');

% %Add geometry information (NodeSet)
% domNode=addGeometryLevel_FEB(domNode,FEB_struct);

%Add boundary condition information
domNode=addBoundaryLevel_FEB(domNode,FEB_struct);

%%
% View example XML string
XML_str = xmlwrite(domNode);
disp(XML_str);

%%

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
