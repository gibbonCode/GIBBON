function [febXML,nodeStruct,elementCell]=import_FEB(febFileName)

% [febXML,nodeStruct,elementCell]=import_FEB(febFileName)
% ------------------------------------------------------------------------
%
% This function imports FEBio XML type input files and generates the
% following outputs: 
%
% febXML : A MATLAB XML object for the input file
% 
% nodeStruct : A structure array containing the following fields: 
%       nodeStruct.N :An nx3 array of nodal coordinates
%       nodeStruct.N_ind :An nx1 array of the node indices (numbers)
%
% elementCell : A cell array containing an entry for each element type
% containing a structure array of the form (for example):
%       E_type: 'quad4'  Specifying element type as the FEBio string
%       E: [4536x4 double] An array of the nodal connectivity
%       E_ind: [4536x1 double] The element indices
%       E_mat: [4536x1 double] The element material indices
%       
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 23/08/2012
%------------------------------------------------------------------------

%%
disp('--- import_FEB ---');

%% IMPORTING .FEB XML
disp('IMPORTING .FEB XML');
febXML = xmlread(febFileName);

%% RETRIEVING NODAL DATA
disp('RETRIEVING NODAL DATA');

GEO_FEB_XML = febXML.getElementsByTagName('Geometry');
NODES_FEB_XML = GEO_FEB_XML.item(0).getElementsByTagName('Nodes');
node_FEB_XML = NODES_FEB_XML.item(0).getElementsByTagName('node');
no_nodes=node_FEB_XML.getLength;

%Allocating memory for nodeStruct
nodeStruct.N=zeros(no_nodes,3);
nodeStruct.N_ind=zeros(no_nodes,1);
for q=0:1:no_nodes-1
    nodeStruct.N_ind(q+1,:)=str2double(node_FEB_XML.item(q).getAttribute('id').toCharArray()');
    nodeStruct.N(q+1,:)=sscanf(node_FEB_XML.item(q).getFirstChild.getData.toCharArray()','%f,%f,%f');
end
disp(['---> Imported ',num2str(no_nodes),' nodes']);

%% RETRIEVING ELEMENT DATA
disp('RETRIEVING ELEMENT DATA');

ELEMENTS_FEB_XML = GEO_FEB_XML.item(0).getElementsByTagName('Elements');
elementTypes={'tri3','quad4','tet4','penta6','hex8'}; %Element type strings, these need to match targets in FEB file, adjust for other types
numberElementNodes=[3,4,4,6,8]; %Number of nodes for each element type

%Getting the number of elements for each type
numberElementEntries=zeros(1,numel(elementTypes));
for q=1:1:numel(elementTypes);
    numberElementEntries(q)=ELEMENTS_FEB_XML.item(0).getElementsByTagName(elementTypes{q}).getLength;
end
indexElementTypes=find(numberElementEntries>0);

%Getting element data
elementType=1;
elementCell=cell(1,numel(indexElementTypes));
for q=indexElementTypes;    
   
    textScanFormat=repmat('%d,',1,numberElementNodes(q)); textScanFormat=textScanFormat(1:end-1); %Create text scan format    
    
    element_FEB_XML=ELEMENTS_FEB_XML.item(0).getElementsByTagName(elementTypes{q});
    
    %Allocating memory for element structure array
    elementStruct.E_type=elementTypes{q}; %Element type
    elementStruct.E=zeros(numberElementEntries(q),str2double(elementTypes{q}(end))); %Element connectivity
    elementStruct.E_ind=zeros(numberElementEntries(q),1); %Element indices
    elementStruct.E_mat=zeros(numberElementEntries(q),1); %Element material indices
    for q2=0:1:numberElementEntries(q)-1 %for each element type
        elementStruct.E(q2+1,:)=sscanf(element_FEB_XML.item(q2).getFirstChild.getData.toCharArray()',textScanFormat); %Element connectivity
        elementStruct.E_ind(q2+1)=str2double(element_FEB_XML.item(q2).getAttribute('id').toCharArray()'); %Element indices
        elementStruct.E_mat(q2+1)=str2double(element_FEB_XML.item(q2).getAttribute('mat').toCharArray()'); %Element material indices
    end    
    elementCell{elementType}=elementStruct; %Add structure to cell array
    elementType=elementType+1; %increase elementTypeIndex
    disp(['---> Imported ',num2str(numberElementEntries(q)),' ',elementTypes{q},' elements']);
end
disp('DONE!');

end

