function [elementStruct,nodeStruct]=import_INP(fileName,numberNodesElement,logicRenumberOption)

% [elementStruct,nodeStruct]=import_INP(fileName,numberNodesElement,logicRenumberOption)
% ------------------------------------------------------------------------
%
% This function imports ABAQUS INP input files and generates the following
% outputs:  
% 
% nodeStruct : A structure array containing the following fields: 
%       nodeStruct.N :An nx3 array of nodal coordinates
%       nodeStruct.N_ind :An nx1 array of the node indices (numbers)
%
% elementStruct : A a structure array of the form (for example):
%       E_type: '*ELEMENT, TYPE=STRI3, ELSET=PART-DEFAULT_1_EB1'  
%               Specifying element type as the ABAQUS string
%       E: [4536x3 double] An mxl array of the nodal connectivity
%       E_ind: [4536x1 double] An mx1 array for the element indices
%
% The users needs to specify the number of nodes in the element using the input:
% numberNodesElement      
% 
% If the input logicRenumberOption is 1 then the nodes are renumbered from
% 1:1:n such that no gaps exist e.g. 2 3 4 7 8 becomes 1 2 3 4 5
%
% The code is developed for importing of single geometry sets only. Model
% boundary conditions and material specifications etc are not imported. The
% case may also not work for multi-part models. 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 23/08/2012
%------------------------------------------------------------------------

%%
disp('--- import_INP ---');

%% IMPORTING .INP FILE INTO CELL ARRAY
disp('IMPORTING .INP FILE INTO CELL ARRAY');
[T]=txtfile2cell(fileName);

%% Finding node and element fields

targets={'*NODE','*ELEMENT','**'}; %Text form targets

%Start scanning
lineCount=1;
lineIndexTarget=zeros(size(targets));
i_target=1;
while 1
    l=T{lineCount};
    target=targets{i_target};
    if (strfind(l,target))
        lineIndexTarget(i_target)=lineCount;
        i_target=i_target+1;
        disp(['---> Found target:',target,' on line ',num2str(lineCount)]);
    end    
    lineCount=lineCount+1;
    if nnz(lineIndexTarget)==(numel(targets))
        break
    end
end
nodeTextField=T(lineIndexTarget(1)+1:lineIndexTarget(2)-3); %NODE TEXT FIELD
elementTextField=T(lineIndexTarget(2)+1:lineIndexTarget(3)-1); %ELEMENT TEXT FIELD

%% CONVERTING TEXT DATA TO MATLAB ARRAYS

disp('CONVERTING TEXT FIELDS TO MATLAB ARRAYS');

%Element data
elementTextScanFormat=repmat('%f,',1,numberNodesElement+1); %+1 due to element index number followed by nodal connectivity
elementTextScanFormat=elementTextScanFormat(1:end-1); %Trim off trailing comma
elementDataField=cell2mat(cellfun(@(x) cell2mat(textscan(x,elementTextScanFormat)),elementTextField,'UniformOutput',0));
E=elementDataField(:,2:end); %The element connectivity matrix
elementIndex=elementDataField(:,1:end); %The element connectivity matrix
elementStruct.E=E;
elementStruct.E_ind=elementIndex(:);
elementStruct.E_type=T{lineIndexTarget(2)};
disp('---> Created elementStruct');
        
%Node data
nodeTextScanFormat=repmat('%f,',1,4); %4 due to node index number followed by x,y,z
nodeTextScanFormat=nodeTextScanFormat(1:end-1); %Trim off trailing comma
nodeDataField=cell2mat(cellfun(@(x) cell2mat(textscan(x,nodeTextScanFormat)),nodeTextField,'UniformOutput',0));
N=nodeDataField(:,2:end); %The nodal coordinates
nodeIndex=nodeDataField(:,1);
nodeStruct.N=N;
nodeStruct.N_ind=nodeIndex(:);
disp('---> Created nodeStruct');

%%

E=elementStruct.E;
N=nodeStruct.N;

%% REMOVING JUMPS IN NODE NUMBERING IF REQUESTED

if logicRenumberOption==1
    nodeIndexNew=1:size(N,1); %The new node indices
    Nn=NaN(max(nodeIndex),3);
    Nn(nodeIndex,:)=N; %contains nan gaps due to jumps
    logicNotUsed=any(isnan(Nn),2); %finds the nan row locations
    nodeIndexSemiFixed=NaN(size(N,1),1); %initialise full of nans
    nodeIndexSemiFixed(~logicNotUsed)=nodeIndexNew; %in valid locations replace by new (other stay nan here but are not used)
    E=nodeIndexSemiFixed(E); %Fixing index numbers    
    elementStruct.E=E; %Overwrite elements
    nodeStruct.N_ind=nodeIndexNew(:); %Overwrite node index arrays  
end
disp('DONE!');
end

