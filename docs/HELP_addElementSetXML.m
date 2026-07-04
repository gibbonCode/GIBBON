%% addElementSetXML
% Below is a demonstration of the features of the |addElementSetXML| function

%%
clear; close all; clc;

%% Syntax
% |[domNode]=addElementSetXML(domNode,rootNode,parseStruct);|

%% Description 
% This function adds an element set defined by input structure parseStruct
% to the XML defined by the XML object domNode.   
% The element value can be a character/string or numerical data.
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
%   <OtherElement>5</OtherElement>
%   <SomeOtherElement>3.1415927e+00</SomeOtherElement>
%   <Material>
%        <material id="1" type="Ogden">
%            <c>2</c>
%        </material>
%   </Material>
% </febio_spec>

%%
% Initialize XML object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec'); 
rootNode=domNode.getElementsByTagName('febio_spec').item(0);

%% Add elements

% Add first element level
parseStruct.OtherElement=5;
parseStruct.SomeOtherElement=pi;
parseStruct.Material=[];
[domNode]=addElementSetXML(domNode,rootNode,parseStruct);

% Add child element of one of the first elements
elementNode=domNode.getElementsByTagName('Material').item(0);
parseStructSub.material=[];
[domNode]=addElementSetXML(domNode,elementNode,parseStructSub);

% Add child element of the child element of one of the firsts elements
subElementNode=domNode.getElementsByTagName('material').item(0);
parseStructSubSub.c=2;
[domNode]=addElementSetXML(domNode,subElementNode,parseStructSubSub);

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
