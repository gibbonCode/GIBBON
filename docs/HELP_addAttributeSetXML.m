%% addAttributeSetXML
% Below is a demonstration of the features of the |addAttributeSetXML| function

%%
clear; close all; clc;

%% Syntax
% |[domNode]=addAttributeSetXML(domNode,elementNode,parseStruct);|

%% Description 
% This function adds a set of attributes defined by the structure
% |parseStruct| the XML defined by the XML object domNode. 
% Attributes are defined in structure arrays. The structure fieldname
% defines the attribute name and the structure value is the attribute
% value. The attribute value can be a character/string or numerical data.
% Vector valued data is transformed to comma separated values e.g. 1:3 will
% become 1, 2, 3. Integer values will remain integer type. For non-integer
% numberical data exponential notation is used e.g. pi becomes
% 3.1415927e+00. 

%% Examples 
%

%% Defining an XML with elements and attributes
% The below example codes this XML section: 
% 
% <febio_spec>
%   <Material>
%       <material id="1" type="Ogden">
%           <c>2</c>
%       </material>
%   </Material>
% </febio_spec>

%%
% Initialize XML object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec'); 
rootNode=domNode.getElementsByTagName('febio_spec').item(0);

%% Add elements
% Main element
elementName='Material';
elementValue=[];
[domNode,elementNode,~]=addElementXML(domNode,rootNode,elementName,elementValue);

%A child element of the main element
elementName='material';
elementValue=[];
[domNode,subElementNode,~]=addElementXML(domNode,elementNode,elementName,elementValue);

%A child element of the child element
elementName='c';
elementValue=2;
[domNode,subSubElementNode,~]=addElementXML(domNode,subElementNode,elementName,elementValue);

%% Defining element attributes

attributeStruct.id=1; %id attribute
attributeStruct.type='Ogden'; %Type attribute
[domNode]=addAttributeSetXML(domNode,subElementNode,attributeStruct);

%% 
% View XML
xmlView(domNode);

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
