%% addGlobalsLevel_FEB
% Below is a demonstration of the features of the |addGlobalsLevel_FEB| function

%%
clear; close all; clc;

%% Syntax
% |[domNode]=addGlobalsLevel_FEB(domNode,FEB_struct);|

%% Description
% This function adds the globals information to the input XML
% data (domNode) based on the input febio structure (FEB_struct).

%% Examples

%% Example: Default globals
% 

%Initialize docNode object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec');

%Add boundary condition information
domNode=addGlobalsLevel_FEB(domNode,[]);

%%
%  View example XML string
XML_str = xmlwrite(domNode);
disp(XML_str);

%% Example: Manually defined globals
% 

FEB_struct.Globals.Constants.Names={'T','R','Fc'};
FEB_struct.Globals.Constants.Entries={0,0,0};

%Initialize docNode object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec');

%Add boundary condition information
domNode=addGlobalsLevel_FEB(domNode,FEB_struct);

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
