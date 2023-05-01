function [domNode]=addAttributeSetXML(domNode,elementNode,parseStruct)

% function [domNode]=addAttributeSetXML(domNode,elementNode,parseStruct)
% -----------------------------------------------------------------------
% This function adds a set of attributes defined by the structure
% |parseStruct| the XML defined by the XML object domNode. 
% Attributes are defined in structure arrays. The structure fieldname
% defines the attribute name and the structure value is the attribute
% value. The attribute value can be a character/string or numerical data.
% Vector valued data is transformed to comma separated values e.g. 1:3 will
% become 1, 2, 3. Integer values will remain integer type. For non-integer
% numberical data exponential notation is used e.g. pi becomes
% 3.1415927e+00. 
%
% 
% -----------------------------------------------------------------------

%%

fieldNameSet = fieldnames(parseStruct);
for q_field=1:1:numel(fieldNameSet)
    attrName=fieldNameSet{q_field};
    attrValue=parseStruct.(attrName);
    [domNode]=addAttributeXML(domNode,elementNode,attrName,attrValue);
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
