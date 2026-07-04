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
% Change log:
% 23/08/2012 Created
% 2014/10/10 Updated for febio_spec 2.0
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
[febio_spec]=getFebioSpecVersion(febXML);

%% RETRIEVING NODAL DATA
[nodeStruct]=get_FEB_nodes(febXML,febio_spec);

%% RETRIEVING ELEMENT DATA

[elementCell]=get_FEB_elements(febXML,febio_spec);

%%

disp('DONE!');

end

%%

function [nodeStruct]=get_FEB_nodes(febXML,febio_spec)

disp('---> Getting nodes');
switch febio_spec
    case '2.0'
        MESH_FEB_XML = febXML.getElementsByTagName('Geometry');
    case '2.5'
        MESH_FEB_XML = febXML.getElementsByTagName('Geometry');
    case '3.0'
        MESH_FEB_XML = febXML.getElementsByTagName('Mesh');
    otherwise %Assume 3.0
        MESH_FEB_XML = febXML.getElementsByTagName('Mesh');
end

NODES_FEB_XML = MESH_FEB_XML.item(0).getElementsByTagName('Nodes');
node_FEB_XML = NODES_FEB_XML.item(0).getElementsByTagName('node');
no_nodes=node_FEB_XML.getLength;

%Allocating memory for nodeStruct
nodeStruct.N=zeros(no_nodes,3);
nodeStruct.N_ind=zeros(no_nodes,1);
for q=0:1:no_nodes-1
    nodeStruct.N_ind(q+1,:)=str2double(node_FEB_XML.item(q).getAttribute('id').toCharArray()');
    nodeStruct.N(q+1,:)=sscanf(node_FEB_XML.item(q).getFirstChild.getData.toCharArray()','%f,%f,%f');
end
disp(['-----> Imported ',num2str(no_nodes),' nodes']);

end

%%

function [elementCell]=get_FEB_elements(febXML,febio_spec)

disp('---> Getting elements');
switch febio_spec
    case {'2.0','2.5'}
        MESH_FEB_XML = febXML.getElementsByTagName('Geometry');
    case '3.0'
        MESH_FEB_XML = febXML.getElementsByTagName('Mesh');
    otherwise %Assume 3.0
        MESH_FEB_XML = febXML.getElementsByTagName('Mesh');
end
Elements_FEB_XML = MESH_FEB_XML.item(0).getElementsByTagName('Elements');


switch febio_spec
    case {'2.0','2.5'}
        numElementsSets=Elements_FEB_XML.getLength; %Number of Elements sets
        elementCell=cell(numElementsSets,1);
        for q=1:1:numElementsSets
            Elements_Set_FEB_XML = Elements_FEB_XML.item(q-1); %The current Elements set
            
            elementStruct.E_type=Elements_Set_FEB_XML.getAttribute('type').toCharArray()';
            elementStruct.E_mat=str2double(Elements_Set_FEB_XML.getAttribute('mat'));
            
            elem_FEB_XML = Elements_Set_FEB_XML.getElementsByTagName('elem'); %The current elem set
            numElem=elem_FEB_XML.getLength; %Number of elem entries in the current Elements set
            numNodesElement=numel(sscanf(elem_FEB_XML.item(0).getFirstChild.getData.toCharArray()','%d,')); %Get num. nodes per element using first
            
            elementStruct.E=zeros(numElem,numNodesElement);
            elementStruct.E_ind=zeros(numElem,1);
            for qe=1:1:numElem
                elementStruct.E(qe,:)=sscanf(elem_FEB_XML.item(qe-1).getFirstChild.getData.toCharArray()','%d,'); %Element nodes
                elementStruct.E_ind(qe)=sscanf(elem_FEB_XML.item(qe-1).getAttribute('id').toCharArray()','%d'); %Element nodes
            end
            elementCell{q}=elementStruct;
            disp(['-----> Imported ',num2str(numElem),' ',elementStruct.E_type,' elements']);
        end
    otherwise %Assume 3.0 
        numElementsSets=Elements_FEB_XML.getLength; %Number of Elements sets
        
        elementCell=cell(numElementsSets,1);
        for q=1:1:numElementsSets
            Elements_Set_FEB_XML = Elements_FEB_XML.item(q-1); %The current Elements set
            
            elementStruct.E_type=Elements_Set_FEB_XML.getAttribute('type').toCharArray()';            
            elementStruct.E_part=Elements_Set_FEB_XML.getAttribute('name').toCharArray()';
            
            elem_FEB_XML = Elements_Set_FEB_XML.getElementsByTagName('elem'); %The current elem set
            numElem=elem_FEB_XML.getLength; %Number of elem entries in the current Elements set
            numNodesElement=numel(sscanf(elem_FEB_XML.item(0).getFirstChild.getData.toCharArray()','%d,')); %Get num. nodes per element using first
            
            elementStruct.E=zeros(numElem,numNodesElement);
            elementStruct.E_ind=zeros(numElem,1);
            for qe=1:1:numElem
                elementStruct.E(qe,:)=sscanf(elem_FEB_XML.item(qe-1).getFirstChild.getData.toCharArray()','%d,'); %Element nodes
                elementStruct.E_ind(qe)=sscanf(elem_FEB_XML.item(qe-1).getAttribute('id').toCharArray()','%d'); %Element nodes
            end            
            elementCell{q}=elementStruct;
            disp(['-----> Imported ',num2str(numElem),' ',elementStruct.E_type,' elements']);
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
