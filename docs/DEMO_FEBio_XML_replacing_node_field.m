%% DEMO_FEBio_XML_replacing_node_field
% Below is a demonstration for:
% 1) Creating basic XML data for a .feb file
% 2) Replacing field in the created XML data by first removing and
% recreating the field.

%%
clear; close all; clc;

%% Creating a feb structure with a geometry section and nodes

%Geometry section
FEB_struct.Geometry.Nodes=rand(8,3); %A set of 8 nodes

%% Using XML coding to build the FEB file

%Initialize docNode object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec');

%Add geometry level to XML file
domNode=addGeometryLevel_FEB(domNode,FEB_struct);

%View the XML data
% xmlView(domNode);
XML_str = gxmlwrite(domNode);
disp(XML_str);

%% Replacing the node field 
% Remove the old node field
geometryNode = domNode.getElementsByTagName('Geometry').item(0); %Get parent node
nodesNode = geometryNode.getElementsByTagName('Nodes').item(0); %Get child node
geometryNode.removeChild(nodesNode); %Remove child node

%%
% Create some new nodes in the struct
FEB_struct.Geometry.Nodes=rand(12,3);

%%
% Add the new nodes to the XML description, see also |febStruct2febFile|

%Create new Nodes child
geometryNode = domNode.getElementsByTagName('Geometry').item(0); %Get parent node
parent_node = domNode.createElement('Nodes'); %Create node node
parent_node = geometryNode.appendChild(parent_node);

n_steps=size(FEB_struct.Geometry.Nodes,1);
for q_n=1:1:n_steps
    node_node = domNode.createElement('node'); %create node entry
    node_node = parent_node.appendChild(node_node); %add node entry
    attr = domNode.createAttribute('id'); %Create id attribute
    attr.setNodeValue(sprintf('%u',q_n)); %Set id text
    node_node.setAttributeNode(attr); %Add id attribute
    node_node.appendChild(domNode.createTextNode(sprintf('%6.7e, %6.7e, %6.7e',FEB_struct.Geometry.Nodes(q_n,:)))); %append data text child
end

%%
%View the new XML data
% xmlView(domNode);
XML_str = gxmlwrite(domNode);
disp(XML_str);

%% Writing the XML data to a .feb file
% Exporting the XML data to a .feb file can be done using |write_XML_no_extra_lines|

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
