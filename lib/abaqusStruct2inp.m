function [varargout]=abaqusStruct2inp(varargin)

% function [T]=abaqusStruct2inp(parseStruct,fileName,optionStruct)
% ------------------------------------------------------------------------
%
% This function provides the basis for coding Abaqus INP input file
% content.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2018/09/06 Created
% 2018/12/04 Added field order entries
% 2020/05/20 Improved performance to fix wrapping. Now uses reshape
% 2020/05/20 Fixed bug in wrapping of non-numeric data
%------------------------------------------------------------------------

%% Parse input

defaultOptionStruct.attributeKeyword='ATTR';
defaultOptionStruct.valueKeyword='VAL';
defaultOptionStruct.commentKeyword='COMMENT';
defaultOptionStruct.customKeyword='CSTM';
defaultOptionStruct.fieldOrder={'heading','preprint','part','node','element',...
                                'surface','distribution','orientation','section',...
                                'assembly','distribution_table','material','Depvar',...
                                'User_Material','contact_pair','surface_interaction',...
                                'step','static','boundary','controls','restart',...
                                'output','element_output','node_output','contact_output'};
defaultOptionStruct.addEnd={'part','assembly','instance','step'};
defaultOptionStruct.addLines=1;

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
    fileName=fullfile(savePath,'temp.inp');
end

file_id = fopen(fileName,'w'); %Open file for writing

%% Reorder fields
[parseStruct]=reorderStruct(parseStruct,optionStruct.fieldOrder);

%% Add content to INP
[file_id]=abaqusStruct2inpStep(file_id,parseStruct,optionStruct);
fclose(file_id); %Close file

%% Collect output
if nargout>0
    T=txtfile2cell(fileName);
    varargout{1}=T;
end

end


function [file_id]=abaqusStruct2inpStep(file_id,parseStruct,optionStruct)

%%
%Get keywords
attributeKeyword=optionStruct.attributeKeyword;
valueKeyword=optionStruct.valueKeyword;
commentKeyword=optionStruct.commentKeyword;
customKeyword=optionStruct.customKeyword;

%Get field names occuring in the input structure
fieldNameSet = fieldnames(parseStruct);

for q_field=1:1:numel(fieldNameSet) %Loop for all field names
    
    loopedOverCell=0;
    
    %Get current field names and values
    currentFieldName=fieldNameSet{q_field}; %attribute name or element name
    currentFieldValue=parseStruct.(fieldNameSet{q_field}); %attribute value or element value
    
    %Initialize field line
    fieldLine=['* ',currentFieldName];
    
    %Replace __ with a hyphen
    fieldLine=regexprep(fieldLine,'__','-');

    %Replace _ with a space    
    fieldLine=strrep(fieldLine,'_',' ');    
    
    if ~isempty(currentFieldValue)
        if ~isstruct(currentFieldValue) %Assume it is value data
            if iscell(currentFieldValue)
                if any(cellfun(@isstruct,currentFieldValue)) %Cell containing structures
                    for q_cell=1:1:numel(currentFieldValue)
                        currentFieldValueSub=currentFieldValue{q_cell};
                        if isa(currentFieldValueSub,'struct')
                            if numel(fieldnames(currentFieldValueSub))~=0 %Not an empty structure
                                optionStructSub=optionStruct;
                                optionStructSub.addLines=0;
                                tempStruct.(currentFieldName)=currentFieldValueSub;
                                [file_id]=abaqusStruct2inpStep(file_id,tempStruct,optionStructSub);
                                clear tempStruct;
                                loopedOverCell=1;
                            end
                        end
                    end
                elseif any(cellfun(@iscell,currentFieldValue)) %Cell containing cells
                    for q_cell=1:1:numel(currentFieldValue)
                        currentFieldValueSub=currentFieldValue{q_cell};
                        optionStructSub=optionStruct;
                        optionStructSub.addLines=0;                        
                        tempStruct.(currentFieldName)=currentFieldValueSub;
                        [file_id]=abaqusStruct2inpStep(file_id,tempStruct,optionStructSub);
                        clear tempStruct;
                    end
                else
                    %Write current field line
                    fprintf(file_id,'%s \n',fieldLine);
                    if numel(currentFieldValue)>1
                        [file_id]=writeValueData(file_id,{currentFieldValue});
                    else
                        [file_id]=writeValueData(file_id,currentFieldValue);
                    end
                end
            else
                %Write current field line
                fprintf(file_id,'%s \n',fieldLine);
                [file_id]=writeValueData(file_id,currentFieldValue);
            end
        else %Assume it is attribute, comment, and value data
            
            %Get current attribute set if any
            if isfield(currentFieldValue,attributeKeyword)
                attributeStruct=currentFieldValue.(attributeKeyword);
                attributeNameSet=fieldnames(attributeStruct);
                for q_attr=1:1:numel(attributeNameSet)
                    attributeNameNow=attributeNameSet{q_attr};
                    attributeValueNow=vec2strIntDouble(attributeStruct.(attributeNameNow),'%6.7e');

                    %Replace __ with a hyphen
                    attributeNameNow=regexprep(attributeNameNow,'__','-');

                    %Replace _ with a space
                    attributeNameNow=strrep(attributeNameNow,'_',' ');

                    if ~isempty(attributeValueNow)
                        fieldLine=[fieldLine,', ',attributeNameNow,'=',attributeValueNow];
                    else %empty attribute field does not use equals sign
                        fieldLine=[fieldLine,', ',attributeNameNow];
                    end
                end
                %Remove attribute field
                currentFieldValue = rmfield(currentFieldValue,attributeKeyword);
            end
            
            %Write current field line with attribute text if any
            fprintf(file_id,'%s \n',fieldLine);
            
            %Add current comments
            if isfield(currentFieldValue,commentKeyword)
                commentData=currentFieldValue.(commentKeyword);
                if iscell(commentData)
                    %Write multiple comment lines
                    for q_comment=1:1:numel(commentData)
                        fprintf(file_id,'%s \n',['** ',commentData{q_comment}]);
                    end
                else
                    %Write comment line
                    fprintf(file_id,'%s \n',['** ',commentData]);
                end
                %Remove comment field
                currentFieldValue = rmfield(currentFieldValue,commentKeyword);
            end

            %Add current empty
            if isfield(currentFieldValue,customKeyword)
                emptyData=currentFieldValue.(customKeyword);
                if iscell(emptyData)
                    %Write multiple empty lines
                    for q_empty=1:1:numel(emptyData)
                        fprintf(file_id,'%s \n',emptyData{q_empty});
                    end
                else
                    %Write empty line
                    fprintf(file_id,'%s \n',emptyData);
                end
                %Remove empty field
                currentFieldValue = rmfield(currentFieldValue,customKeyword);
            end
            
            %Add current values
            if isfield(currentFieldValue,valueKeyword)
                [file_id]=writeValueData(file_id,{currentFieldValue.(valueKeyword)});
                
                %Remove value field
                currentFieldValue = rmfield(currentFieldValue,valueKeyword);
            end
            
            %Recursively parse nested remainder of structure
            if numel(fieldnames(currentFieldValue))~=0 %Not an empty structure
                optionStructSub=optionStruct;
                optionStructSub.addLines=0;
                [file_id]=abaqusStruct2inpStep(file_id,currentFieldValue,optionStructSub);
            end
        end
    end
    if any(strcmpi(optionStruct.addEnd,currentFieldName)) && loopedOverCell==0
        fprintf(file_id,'%s \n',['* End ',currentFieldName]); %Also write end statement for master field
    end
    
    if optionStruct.addLines==1
        fprintf(file_id,'%s \n',['** ',repmat('-',1,100)]);
    end
end

end

%%
function [B]=reorderStruct(A,fieldOrder)

%Start with fieldOrder fields
fieldNameSet=fieldnames(A);
for q=1:1:numel(fieldOrder)
    fieldNameNow=fieldOrder{q};
    logicFind=strcmpi(fieldNameSet,fieldNameNow);
    if any(logicFind)
        B.(fieldNameSet{logicFind})=A.(fieldNameSet{logicFind});
        A=rmfield(A,fieldNameSet{logicFind});
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
function [file_id]=writeValueData(file_id,valueData)
if iscell(valueData)
    %Write multiple value entries
    for q_value=1:1:numel(valueData)
        
        if iscell(valueData{q_value})
            C=valueData{q_value};
            if iscolumn(C)
                for qr=1:1:numel(C)
                   c=C{qr} ;
                    for q_value_col=1:1:numel(c)
                        if q_value_col==1
                            t = vec2strIntDouble(c{q_value_col},'%6.7e');
                            t=strsplit(t,'\n')';
                        else
                            tn = vec2strIntDouble(c{q_value_col},'%6.7e');
                            TN=strsplit(tn,'\n')';
                            D=repmat({', '},size(t,1),1);
                            t = strcat(t,D,TN);
                        end
                    end
                    for qt=1:1:numel(t)
                        fprintf(file_id,'%s \n',t{qt});
                    end
                end
            else
                for q_value_col=1:1:numel(C)
                    if q_value_col==1
                        t = vec2strIntDouble(C{q_value_col},'%6.7e');
                        t=strsplit(t,'\n')';
                    else
                        tn = vec2strIntDouble(C{q_value_col},'%6.7e');
                        TN=strsplit(tn,'\n')';
                        D=repmat({', '},size(t,1),1);
                        t = strcat(t,D,TN);
                    end
                end
                for qt=1:1:numel(t)
                    fprintf(file_id,'%s \n',t{qt});
                end
            end
            
        else
            t=toTextCheckWrap(valueData{q_value},16); 
            fprintf(file_id,'%s \n',t);
        end
        
    end
else
    %Write value entry
    t=toTextCheckWrap(valueData,16);
    fprintf(file_id,'%s \n',t);
end
end

%%

function [t]=toTextCheckWrap(valueData,wrapLength)

if isnumeric(valueData)
    if isrow(valueData) && numel(valueData)>wrapLength
        %OLD and slow        
        %t=vec2strIntDouble(valueData,'%6.7e');
        %t=strwrap(t,wrapLength,', '); %Wrap to max width of wrapLength entries

        if rem(numel(valueData),wrapLength)>0
            valueDataTemp=nan(1,ceil(numel(valueData)/wrapLength)*wrapLength);
            valueDataTemp(1:numel(valueData))=valueData;
            valueDataTemp=reshape(valueDataTemp,wrapLength,numel(valueDataTemp)/wrapLength)';
            t=vec2strIntDouble(valueDataTemp,'%6.7e');
            t=regexprep(t,',(\s*)NaN',''); %remove nan
        else
            valueData=reshape(valueData,wrapLength,numel(valueData)/wrapLength)';
            t=vec2strIntDouble(valueData,'%6.7e');
        end                
    else
        t=vec2strIntDouble(valueData,'%6.7e');
    end
else
    t=valueData;
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
