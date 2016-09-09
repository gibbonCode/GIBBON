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