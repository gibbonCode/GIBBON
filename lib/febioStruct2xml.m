function [varargout]=febioStruct2xml(varargin)

% function [domNode]=febioStruct2xml(parseStruct,fileName,optionStruct)
% ------------------------------------------------------------------------
%
% This function provides the basis for coding FEBio XML input file content
% and helps to generate the XML input file.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2018/01/15 Created
% 2018/02/09 Found bug in text write mode.
% 2018/02/18 Created several "array parse methods" including based on text
% file exporting providing a significant computational time reduction over
% full XML based parsing. 
% 2018/05/15 Create temp directory if it is does not exist
% To do:
% Gracefully handle empty fields e.g. 0x8 element array
%------------------------------------------------------------------------

%% Parse input

defaultOptionStruct.attributeKeyword='ATTR';
defaultOptionStruct.valueKeyword='VAL';
defaultOptionStruct.arrayLoopKeywords={'node','elem','face','delem','quad4','quad8','tri3','tri6','tri7','line2','line3','point'};
defaultOptionStruct.fieldOrder={'Module','Control','Globals','Material',...
    'Geometry','MeshData','Initial','Boundary',...
    'Loads','Contact','Constraints','Rigid Connectors',...
    'Discrete','LoadData','Output','Parameters'};
defaultOptionStruct.arrayParseMethod=1;

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

defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

%Create temp folder if it does not exist
if ~exist(savePath,'file')
    mkdir(savePath);
end

%Create output file name if not provided
if isempty(fileName)
    fileName=fullfile(savePath,'temp.xml');
end

%% Initialize XML object
domNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec'); %Create the overall febio_spec field

%% Add comment
commentString = ['Created using GIBBON, ',datestr(now)];
commentNode = domNode.createComment(commentString);
rootNode=domNode.getElementsByTagName('febio_spec').item(0);
rootNode.appendChild(commentNode);

%% Reorder fields
[parseStruct]=reorderStruct(parseStruct,optionStruct.fieldOrder);

%% Add content to XML
[domNode]=febioStruct2xmlStep(domNode,rootNode,parseStruct,optionStruct,fileName);

%% Export XML file
xmlwrite_xerces(fileName,domNode); %Custom XML write function

%Import to cut out bad element names from text write mode
T=txtfile2cell(fileName);

logicKeep=~gcontains(T,'xxxxxxxxxxx_CUT');
cell2txtfile(fileName,T(logicKeep),0,0);

%% Collect output
if nargout>0
    varargout{1}=domNode;
end

end

function [domNode,rootNode]=febioStruct2xmlStep(domNode,rootNode,parseStruct,optionStruct,fileName)

arrayParseMethod=optionStruct.arrayParseMethod;

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
                    switch arrayParseMethod
                        case 0
                            for q_array=1:1:numArray %Loop over elements in array
                                %Add the element and its value
                                if iscell(arrayDataSet)
                                    [domNode,elementNode,rootNode]=addElementXML(domNode,rootNode,currentFieldName,arrayDataSet{q_array});
                                else
                                    [domNode,elementNode,rootNode]=addElementXML(domNode,rootNode,currentFieldName,arrayDataSet(q_array,:));
                                end
                                if ~isempty(fieldNameSetSub) %If attributes exist
                                    for q_fieldSub=1:1:numel(fieldNameSetSub) %loop over all attributes
                                        currentFieldNameSub=fieldNameSetSub{q_fieldSub}; %Get current attribute name
                                        %Create structure for attribute creation
                                        if iscell(attributeStruct.(currentFieldNameSub))
                                            attributeStructSub.(currentFieldNameSub)=attributeStruct.(currentFieldNameSub){q_array};
                                        else
                                            attributeStructSub.(currentFieldNameSub)=attributeStruct.(currentFieldNameSub)(q_array,:);
                                        end
                                        [domNode]=addAttributeSetXML(domNode,elementNode,attributeStructSub); %Add attribute
                                    end
                                end
                            end
                        otherwise
                            [domNode,rootNode]=textModeElementAttributeParse(domNode,rootNode,currentFieldName,attributeStruct,arrayDataSet,numArray,fieldNameSetSub,fileName,arrayParseMethod);
                    end
                    
                elseif ~isempty(fieldNameSetSub) %If this is not empty there are attributes to add despite empty array
                    
                    %Get length from first
                    currentFieldNameSub=fieldNameSetSub{1};
                    numArray=numel(attributeStruct.(currentFieldNameSub));
                    
                    switch arrayParseMethod
                        case 0
                            for q_array=1:1:numArray %Loop over all entries
                                [domNode,elementNode,rootNode]=addElementXML(domNode,rootNode,currentFieldName,[]); %Create element without value
                                for q_fieldSub=1:1:numel(fieldNameSetSub) %Loop over and add all attributes
                                    currentFieldNameSub=fieldNameSetSub{q_fieldSub};
                                    %Create structure for attribute creation
                                    if iscell(attributeStruct.(currentFieldNameSub))
                                        attributeStructSub.(currentFieldNameSub)=attributeStruct.(currentFieldNameSub){q_array};
                                    else
                                        attributeStructSub.(currentFieldNameSub)=attributeStruct.(currentFieldNameSub)(q_array,:);
                                    end
                                    [domNode]=addAttributeSetXML(domNode,elementNode,attributeStructSub); %Add attribute
                                end
                            end
                        otherwise
                            [domNode,rootNode]=textModeElementAttributeParse(domNode,rootNode,currentFieldName,attributeStruct,arrayDataSet,numArray,fieldNameSetSub,fileName,arrayParseMethod);
                    end
                elseif ~isempty(parseStructSub) %Despite arrayLoopKeywords trigger no values or attributes are to be added, but other fields might exist
                    [domNode,rootNode]=recursiveElementParse(domNode,rootNode,currentFieldName,currentFieldValue,optionStruct,fileName);
                end
                
            else
                [domNode,rootNode]=recursiveElementParse(domNode,rootNode,currentFieldName,currentFieldValue,optionStruct,fileName);
            end
    end
end

end

function [domNode,rootNode]=recursiveElementParse(domNode,rootNode,currentFieldName,currentFieldValue,optionStruct,fileName)

if isnumeric(currentFieldValue) || ischar(currentFieldValue)
    [domNode,~,rootNode]=addElementXML(domNode,rootNode,currentFieldName,currentFieldValue);
elseif isstruct(currentFieldValue)
    [domNode,elementNode,rootNode]=addElementXML(domNode,rootNode,currentFieldName,[]); %Create the current element
    [domNode,elementNode]=febioStruct2xmlStep(domNode,elementNode,currentFieldValue,optionStruct,fileName); %Add whatever is nested in here
elseif iscell(currentFieldValue)
    for q_cell=1:1:numel(currentFieldValue)
        parseStructSub=currentFieldValue{q_cell};
        [domNode,elementNode,rootNode]=addElementXML(domNode,rootNode,currentFieldName,[]);
        [domNode,elementNode]=febioStruct2xmlStep(domNode,elementNode,parseStructSub,optionStruct,fileName); %Add whatever is nested in here
    end
end

end

%%

function [domNode,rootNode]=textModeElementAttributeParse(domNode,rootNode,currentFieldName,attributeStruct,arrayDataSet,numArray,fieldNameSetSub,fileName,arrayParseMethod)

switch arrayParseMethod
    case 1 %Non-flexible, assumes equal attribute type and size (fastest)
        textFormat=['<',currentFieldName,' '];
        if ~isempty(fieldNameSetSub) %If attributes exist
            D=[];
            for q_fieldSub=1:1:numel(fieldNameSetSub) %loop over all attributes
                attributeName=fieldNameSetSub{q_fieldSub}; %Get current attribute name
                if iscell(attributeStruct.(attributeName))
                    attributeNameValue=attributeStruct.(attributeName);
                else
                    attributeNameValue=attributeStruct.(attributeName);
                end
                D=[D attributeNameValue];
                if all(isrounded(attributeNameValue))
                    t_form=repmat(['%d',', '],1,size(attributeNameValue,2));
                else
                    t_form=repmat(['%6.7e',', '],1,size(attributeNameValue,2));
                end
                t_form=t_form(1:end-2); %Take away last space and comma
                textFormat=[textFormat,' ',attributeName,'="',t_form,'"'];
            end
        else
            D=[];
        end
        
        if all(isrounded(arrayDataSet))
            t_form=repmat(['%d',', '],1,size(arrayDataSet,2));
        else
            t_form=repmat(['%6.7e',', '],1,size(arrayDataSet,2));
        end
        t_form=t_form(1:end-2); %Take away last space and comma
        textFormat=[textFormat,'>',t_form,'</',currentFieldName,'>','\n'];
        D=[D arrayDataSet];
        
        %Convert to string
        t=sprintf(textFormat,D');
        
        T_add={t};      
    case 2 %Flexible, handles non-uniform attributes, uses cellfun to "loop"         
        T_add=cell(numArray,1);
        tquote=repmat({'"'},numArray,1);
        if ~isempty(fieldNameSetSub) %If attributes exist
            for q_fieldSub=1:1:numel(fieldNameSetSub) %loop over all attributes
                attributeName=fieldNameSetSub{q_fieldSub}; %Get current attribute name
                attributeNameValue=attributeStruct.(attributeName);                
                T_attributeValue=var2cellstr(attributeNameValue);          
                T_add=strcat(T_add,repmat({[' ',attributeName,'=']},numArray,1),tquote,T_attributeValue,tquote);
            end
        end
        T_val=var2cellstr(arrayDataSet);
        T_add=strcat(repmat({['<',currentFieldName,' ']},numArray,1),T_add,repmat({'>'},numArray,1),T_val,repmat({['</',currentFieldName,'>']},numArray,1));    
    case 3 %Flexible, handles non-uniform attributes, loop through all entries (slowest)
        T_add=cell(numArray,1);
        for q_array=1:1:numArray %Loop over elements in array
            attributeText=[];
            if ~isempty(fieldNameSetSub) %If attributes exist
                for q_fieldSub=1:1:numel(fieldNameSetSub) %loop over all attributes                    
                    attributeName=fieldNameSetSub{q_fieldSub}; %Get current attribute name
                    if iscell(attributeStruct.(attributeName))
                        attributeNameValue=attributeStruct.(attributeName){q_array};
                    else
                        attributeNameValue=attributeStruct.(attributeName)(q_array,:);
                    end                    
                    [attributeNameValue]=vec2strIntDouble(attributeNameValue,'%6.7e');                    
                    attributeText=[attributeText,' ',attributeName,'="',attributeNameValue,'"'];
                end
            end
            if ~isempty(arrayDataSet)
                [currentValue]=vec2strIntDouble(arrayDataSet(q_array,:),'%6.7e');
            else
                currentValue=[];
            end
            T_add{q_array}=['<',currentFieldName,' ',attributeText,'>',currentValue,'</',currentFieldName,'>'];
        end
end

fieldNameCut=[currentFieldName,'_xxxxxxxxxxx_CUT'];
T_top{1}=['<',fieldNameCut,'>'];
totalTextCell=T_top;
totalTextCell(end+1:end+numel(T_add))=T_add;
totalTextCell(end+1)={['</',fieldNameCut,'>']};

% Export text cell to .feb file
cell2txtfile(fileName,totalTextCell,0,0); %Export to text file

% Reimport XML type
domNode2 = xmlread(fileName); %Reimport docNode

nodeList=domNode2.getElementsByTagName(fieldNameCut);
rootNode2=nodeList.item(0);

firstDocImportedNode = domNode.importNode(rootNode2.getChildNodes,true);
rootNode.appendChild(firstDocImportedNode);

%Delete file
delete(fileName);

end

%%
function [B]=reorderStruct(A,fieldOrder)

%Start with fieldOrder fields
fieldNameSet=fieldnames(A);
for q=1:1:numel(fieldOrder)
    fieldNameNow=fieldOrder{q};
    logicFind=strcmp(fieldNameNow,fieldNameSet);
    if any(logicFind)
        B.(fieldNameNow)=A.(fieldNameNow);
        A=rmfield(A,fieldNameNow);
    end
end

%Add remaining fields
fieldNameSet=fieldnames(A);
for q=1:1:numel(fieldNameSet)
    fieldNameNow=fieldNameSet{q};
    B.(fieldNameNow)=A.(fieldNameNow);
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
% Copyright (C) 2019  Kevin Mattheus Moerman
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
