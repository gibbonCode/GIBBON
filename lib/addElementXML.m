function [varargout]=addElementXML(domNode,rootNode,elementName,elementValue,varargin)


%%
try
    elementNode = domNode.createElement(elementName); %create entry
    elementNode = rootNode.appendChild(elementNode); %add entry
catch
    % Get correct rootNode
    nodeList=domNode.getElementsByTagName(rootNode.getNodeName);
    rootNode=nodeList.item(nodeList.getLength-1);
    
    elementNode = domNode.createElement(elementName); %create entry
    elementNode = rootNode.appendChild(elementNode); %add entry    
end

if ~isempty(elementValue)    
    %Parse optional option struct input
    defaultOptionStruct.formatDouble='%6.7e';
    defaultOptionStruct.formatInteger='%d';
    defaultOptionStruct.dlmChar=',';
    defaultOptionStruct.rowWrapLength=[];

    if nargin==5
        optionStruct=varargin{1};
        [optionStruct]=structComplete(optionStruct,defaultOptionStruct,1); %Complement provided with default if missing or empty
    else
        optionStruct=defaultOptionStruct;
    end

    [elementValue]=mat2strIntDouble(elementValue,optionStruct);
    elementNode.appendChild(domNode.createTextNode(elementValue)); %append data text child
end

%Collect output
varargout{1}=domNode;
varargout{2}=elementNode;
varargout{3}=rootNode;

end

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
