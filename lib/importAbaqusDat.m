function [abaqusData]=importAbaqusDat(varargin)

% function [abaqusData]=importAbaqusDat(fileName,optionStruct)
%-------------------------------------------------------------------------
%
% Form:
%         NODE FOOT-  COOR1          COOR2          COOR3
%         NOTE
%
%         1     -5.0000000E+00 -5.0000000E+00 -5.0000000E+00
%         2     -5.0000000E+00 -1.1555553E+00 -5.0000000E+00
%         ....
%         n     -5.0000000E+00 -1.1555553E+00 -5.0000000E+00
%
%         MAXIMUM          6.533       6.533       2.000
%         AT NODE                13          20          49
%
%         MINIMUM         -5.000      -5.000      -5.000
%         AT NODE                 1           1           111

% Form:
%         ELEMENT  PT FOOT-       E11         E22         E33         E12         E13         E23
%         NOTE
%
%         1   1         0.1427      0.1427     -0.3566      3.0358E-18  1.9120E-17  1.7250E-18
%         1   2         0.1427      0.1427     -0.3566      2.0817E-17  3.9274E-17 -2.5351E-17
%         ....
%         n   8         0.1427      0.1427     -0.3566      6.5052E-17  3.6108E-17 -4.4747E-17
%
%         MAXIMUM                0.1427      0.1427     -0.3566      9.8446E-17  2.0358E-16  2.1077E-16
%         ELEMENT                     6           1           3          24           7          15
%
%         MINIMUM                0.1427      0.1427     -0.3566     -8.8471E-17 -1.6717E-16 -1.7420E-16
%         ELEMENT                     8           3          10           2          12          11
%
%-------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        fileName=varargin{1};
        optionStruct=[];
    case 2
        fileName=varargin{1};
        optionStruct=varargin{2};
end

defaultOptionStruct.nodeData.dataTypeStr='N O D E   O U T P U T';
defaultOptionStruct.nodeData.startTarget='    NODE FOOT-';
defaultOptionStruct.elementData.dataTypeStr='E L E M E N T   O U T P U T';
defaultOptionStruct.elementData.startTarget='    ELEMENT  PT';
defaultOptionStruct.stopTarget=' MAXIMUM';
defaultOptionStruct.getElementData=1;
defaultOptionStruct.getNodeData=1;

[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1); %Complement provided with default if missing or empty

%%

T=txtfile2cell(fileName);

for qStep=1:1:1
    
    incrementLines = find(~cellfun(@isempty,regexp(T,'INCREMENT.*SUMMARY')));
    incrementLines(end+1)=numel(T);
    
    for qIncrement=1:1:numel(incrementLines)-1
        
        %Get increment data
        s=T{incrementLines(qIncrement)};
        %Scans this type of line: '                                INCREMENT     1 SUMMARY'
        abaqusData.STEP(qStep).INCREMENT(qIncrement).id=str2double(cell2mat(regexp(s,'\d+\.?\d*|-\d+\.?\d*|\.?\d+|-\.?\d+','match')));
        
        %Scans this type of line: ' TIME INCREMENT COMPLETED  5.000E-02,  FRACTION OF STEP COMPLETED  5.000E-02'
        s=T{incrementLines(qIncrement)+3};
        d=cellfun(@str2double,regexp(s,'([\+-]?((\d*\.\d+)|\d+))(E[\+-]?\d+)?','match'));
        abaqusData.STEP(qStep).INCREMENT(qIncrement).TIME_INCREMENT_COMPLETED=d(1);
        abaqusData.STEP(qStep).INCREMENT(qIncrement).FRACTION_OF_STEP_COMPLETED=d(2);
        
        %Scans this type of line: '  STEP TIME COMPLETED       5.000E-02,  TOTAL TIME COMPLETED        5.000E-02'
        s=T{incrementLines(qIncrement)+4};
        d=cellfun(@str2double,regexp(s,'([\+-]?((\d*\.\d+)|\d+))(E[\+-]?\d+)?','match'));
        abaqusData.STEP(qStep).INCREMENT(qIncrement).STEP_TIME_COMPLETED=d(1);
        abaqusData.STEP(qStep).INCREMENT(qIncrement).TOTAL_TIME_COMPLETED=d(2);
        
        %Get current data containing part
        T_sub=T(incrementLines(qIncrement):incrementLines(qIncrement+1)-1);
        
        %Get data
        stopTarget=optionStruct.stopTarget;
        if optionStruct.getNodeData==1 %Get node data
            dataTypeStr=optionStruct.nodeData.dataTypeStr;
            startTarget=optionStruct.nodeData.startTarget;
            try
                abaqusData.STEP(qStep).INCREMENT(qIncrement).nodeOutput=getDataField(T_sub,dataTypeStr,startTarget,stopTarget);
            catch
                abaqusData.STEP(qStep).INCREMENT(qIncrement).nodeOutput=NaN;
            end 
        end
        if optionStruct.getElementData==1 %Get element data
            dataTypeStr=optionStruct.elementData.dataTypeStr;
            startTarget=optionStruct.elementData.startTarget;
            try
                abaqusData.STEP(qStep).INCREMENT(qIncrement).elementOutput=getDataField(T_sub,dataTypeStr,startTarget,stopTarget);
            catch
                abaqusData.STEP(qStep).INCREMENT(qIncrement).elementOutput=NaN;
            end
        end                
    end
end
end

%%

function S=getDataField(T_sub,dataTypeStr,startTarget,stopTarget)

dataOutput = any(~cellfun(@isempty,regexp(T_sub,dataTypeStr)));
if dataOutput==1 %Element data found
    dataOutputStartLines = find(~cellfun(@isempty,regexp(T_sub,startTarget))); %Start lines
    dataOutputEndLines = find(~cellfun(@isempty,regexp(T_sub,stopTarget)))-2; %End lines
    dataOutputEndLines=dataOutputEndLines(dataOutputEndLines>dataOutputStartLines(1));
    ignoreSet={'FOOT-'};
    for q_dataOutput=1:1:numel(dataOutputStartLines)
        s=T_sub{dataOutputStartLines(q_dataOutput)}; %Header line
        s = regexprep(s,'\s+',' '); %Remove excessive spaces
        s=s(2:end-1); %Remove initial and final space        
        headerSet=strsplit(s,' '); %Header set in cell        
        %         numHeadersAll=numel(headerSet);
        logicRemove=gcontains(headerSet,ignoreSet);
        headerSet=headerSet(~logicRemove);
        numHeadersKeep=numel(headerSet);
        s=T_sub(dataOutputStartLines(q_dataOutput)+3:dataOutputEndLines(q_dataOutput)); %Data text field        
        try
            D=cell2mat(cellfun(@(x) sscanf(x,'%f')',s,'uniformOutput',0)); %Data array
        catch %Alternative which ignores character components
            D=cell2mat(cellfun(@str2double,regexp(s,'([\+-]?((\d*\.\d+)|\d+))(E[\+-]?\d+)?','match'),'UniformOutput',0)); %Data array
        end
        
        %Create output structure with subfields names after headers
        for q=1:1:numHeadersKeep
            fieldNameNow=headerSet{q} ;
            S(q_dataOutput).data.(fieldNameNow)=D(:,q);
        end
        
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
