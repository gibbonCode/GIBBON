function [varargout]=febioStruct2xml(varargin)

% function [domNode]=febioStruct2xml(parseStruct,fileName,optionStruct)

%% Parse input

defaultOptionStruct.attributeKeyword='ATTR';
defaultOptionStruct.valueKeyword='VAL';
defaultOptionStruct.arrayLoopKeywords={'node','elem','face','delem','quad4','quad8','tri3','tri6','tri7','line2','line3','point'};
defaultOptionStruct.arrayWriteMode='text'; %'text' or 'xml'

switch nargin
    case 1
        parseStruct=varargin{1};
        fileName=[];
        optionStruct=[];
    case 2        
        parseStruct=varargin{1};
        fileName=varargin{2};                
        optionStruct=[];
    case 3
        parseStruct=varargin{1};
        fileName=varargin{2};                
        optionStruct=varargin{3};
    otherwise
        error('Wrong number of intput arguments, between 1 and 3 input arguments supported i.e.: parseStruct,fileName,optionStruct')
end

[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1); %Complement provided with default if missing or empty

%% Initialize XML object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec'); %Create the overall febio_spec field

%% Add comment
commentString = ['Created using GIBBON, ',datestr(now)];
commentNode = domNode.createComment(commentString);
rootNode=domNode.getElementsByTagName('febio_spec').item(0);
rootNode.appendChild(commentNode);

%% Add content to XML
[domNode]=febioStruct2xmlStep(domNode,rootNode,parseStruct,optionStruct);

%% Export XML file
if ~isempty(fileName)
    exportFEB_XML(fileName,domNode);
end

%% Collect output
if nargout>0
   varargout{1}=domNode;
end

end

function [domNode]=febioStruct2xmlStep(domNode,rootNode,parseStruct,optionStruct)

writeMode=1; 

%%
%Get keywords
attributeKeyword=optionStruct.attributeKeyword;
valueKeyword=optionStruct.valueKeyword;
arrayLoopKeywords=optionStruct.arrayLoopKeywords;

%Get field names occuring in the input structure
fieldNameSet = fieldnames(parseStruct);

for q_field=1:1:numel(fieldNameSet) %Loop for all field names
    
    %Get current field names and values
    currentFieldName=fieldNameSet{q_field}; %attribute name or element name
    currentFieldValue=parseStruct.(fieldNameSet{q_field}); %attribute value or element value
    
    %Switch behaviour depending on name, i.e. create attribute, assign
    %value (to current entity), or create an element with a value
    switch currentFieldName
        case attributeKeyword %ATTRIBUTE
            attributeStruct=currentFieldValue; %The attribute value should be a structure e.g. ATTR.attributeName.attributeValue
            if ~isempty(rootNode)
                %Add the current set of attribute for the element at rootNode
                [domNode]=addAttributeSetXML(domNode,rootNode,attributeStruct);
            else
                error('The rootNode is empty, the current attribute appears to lack a parent element');
            end
        case valueKeyword %VALUE
            %Add the current value to rootNode
            [domNode]=addElementValueXML(domNode,rootNode,currentFieldValue);
        otherwise %ELEMENT
            
            %Check if the current element tag correspods to one requiring
            %looping
            logicArray=any(strcmp(currentFieldName,arrayLoopKeywords));
            
            if logicArray %Check for looping across values
                parseStructSub=currentFieldValue; %The current element should be a structure
                if isfield(parseStructSub,valueKeyword) %if the structure contains a fielname which mathes the value keyword
                    arrayDataSet=parseStructSub.(valueKeyword); %Get the array of data
                    parseStructSub = rmfield(parseStructSub,valueKeyword); %Remove the value field from the array
                    numArray=size(arrayDataSet,1); %Measure the number of entries in the array
                else %Some other field
                    arrayDataSet=[];
                    numArray=[];
                end
                
                if isfield(parseStructSub,attributeKeyword) %if the structure contains a fielname which mathes the attribute keyword
                    attributeStruct=parseStructSub.(attributeKeyword);
                    parseStructSub = rmfield(parseStructSub,attributeKeyword); %Remove the value field from the array
                    fieldNameSetSub=fieldnames(attributeStruct);
                else
                    fieldNameSetSub=[];
                end
                
                if ~isempty(arrayDataSet) %if array entries exist
                    switch writeMode
                        case 1
                            for q_array=1:1:numArray %Loop over elements in array
                                %Add the element and its value
                                [domNode,elementNode]=addElementXML(domNode,rootNode,currentFieldName,arrayDataSet(q_array,:));
                                
                                if ~isempty(fieldNameSetSub) %If attributes exist
                                    for q_fieldSub=1:1:numel(fieldNameSetSub) %loop over all attributes
                                        currentFieldNameSub=fieldNameSetSub{q_fieldSub}; %Get current attribute name
                                        attributeStructSub.(currentFieldNameSub)=attributeStruct.(currentFieldNameSub)(q_array,:); %Create structure for attribute creation
                                        [domNode]=addAttributeSetXML(domNode,elementNode,attributeStructSub); %Add attribute
                                    end
                                end
                            end
%                         case 2
%                             T=cell(numArray,1);
%                             for q_array=1:1:numArray %Loop over elements in array
%                                 if ~isempty(fieldNameSetSub) %If attributes exist
%                                     attributeText=[];
%                                     for q_fieldSub=1:1:numel(fieldNameSetSub) %loop over all attributes
%                                         attributeName=fieldNameSetSub{q_fieldSub}; %Get current attribute name
%                                         attributeNameValue=attributeStruct.(currentFieldNameSub)(q_array,:);
%                                         
%                                         attributeText=[attributeText,' ',currentFieldNameSub,'="','" '];
%                                     end
%                                 end
%                                 textLine=['<',currentFieldName,' ',attributeText,'>',]
%                             end
%                             fdsfsa
                    end
                    
                elseif ~isempty(fieldNameSetSub) %If this is not empty there are attributes to add despite empty array
                    %Get length from first
                    currentFieldNameSub=fieldNameSetSub{1};
                    numArray=numel(attributeStruct.(currentFieldNameSub));
                    for q_array=1:1:numArray %Loop over all entries
                        [domNode,elementNode]=addElementXML(domNode,rootNode,currentFieldName,[]); %Create element without value
                        for q_fieldSub=1:1:numel(fieldNameSetSub) %Loop over and add all attributes
                            currentFieldNameSub=fieldNameSetSub{q_fieldSub};
                            attributeStructSub.(currentFieldNameSub)=attributeStruct.(currentFieldNameSub)(q_array,:); %Create structure for attribute creation
                            [domNode]=addAttributeSetXML(domNode,elementNode,attributeStructSub); %Add attribute
                        end
                    end
                elseif ~isempty(parseStructSub) %Despite arrayLoopKeywords trigger no values or attributes are to be added, but other fields might exist
                    [domNode]=recursiveElementParse(domNode,rootNode,currentFieldName,currentFieldValue,optionStruct);
                end
                
            else
                [domNode]=recursiveElementParse(domNode,rootNode,currentFieldName,currentFieldValue,optionStruct);
            end
    end
end

end

function [domNode]=recursiveElementParse(domNode,rootNode,currentFieldName,currentFieldValue,optionStruct)

if isnumeric(currentFieldValue) || ischar(currentFieldValue)
    [domNode,~]=addElementXML(domNode,rootNode,currentFieldName,currentFieldValue);
elseif isstruct(currentFieldValue)
    [domNode,elementNode]=addElementXML(domNode,rootNode,currentFieldName,[]); %Create the current element
    [domNode]=febioStruct2xmlStep(domNode,elementNode,currentFieldValue,optionStruct); %Add whatever is nested in here
elseif iscell(currentFieldValue)
    for q_cell=1:1:numel(currentFieldValue)
        parseStructSub=currentFieldValue{q_cell};
        [domNode,elementNode]=addElementXML(domNode,rootNode,currentFieldName,[]);
        [domNode]=febioStruct2xmlStep(domNode,elementNode,parseStructSub,optionStruct); %Add whatever is nested in here
    end
end

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
% Copyright (C) 2017  Kevin Mattheus Moerman
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