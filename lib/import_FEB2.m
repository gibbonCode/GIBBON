function [T,FEB_XML]=import_FEB2(filename)

%Import .feb file from FEBio
%
%
%23/08/2012, Kevin Mattheus Moerman

%% Loading .feb
[T]=txtfile2cell(filename);
FEB_XML = xmlread(filename);

%% Finding material
MAT_FEB_XML = FEB_XML.getElementsByTagName('Material');
mat_FEB_XML = MAT_FEB_XML.item(0).getElementsByTagName('material');
no_mats=mat_FEB_XML.getLength;
GEO_FEB_XML = FEB_XML.getElementsByTagName('Geometry');

%% Retrieving nodal data

NODES_FEB_XML = GEO_FEB_XML.item(0).getElementsByTagName('Nodes');
node_FEB_XML = NODES_FEB_XML.item(0).getElementsByTagName('node');
no_nodes=node_FEB_XML.getLength;

N=zeros(no_nodes,3);
N_ind=zeros(no_nodes,1);
for i=0:1:no_nodes-1
    N_ind(i+1,:)=str2double(node_FEB_XML.item(i).getAttribute('id').toCharArray()');
    N(i+1,:)=sscanf(node_FEB_XML.item(i).getFirstChild.getData.toCharArray()','%f,%f,%f');
end

%% Retrieving element data

ELEMENTS_FEB_XML = GEO_FEB_XML.item(0).getElementsByTagName('Elements');
element_types={'tri3','quad4','tet4','penta6','hex8'};
text_scan_formats={'%f,%f,%f','%f,%f,%f,%f','%f,%f,%f,%f','%f,%f,%f,%f,%f,%f','%f,%f,%f,%f,%f,%f,%f,%f'};

no_elements=zeros(1,numel(element_types));
for i=1:1:numel(element_types);
    no_elements(i)=ELEMENTS_FEB_XML.item(0).getElementsByTagName(element_types{i}).getLength;
end
IND_element_types=find(no_elements>0);

el_type=1;
E_cell=cell(1,numel(IND_element_types));
for i=IND_element_types;
    element_FEB_XML=ELEMENTS_FEB_XML.item(0).getElementsByTagName(element_types{i});
    
    %Allocating memory
    E_struct.E_type=element_types{i};
    E_struct.E=zeros(no_elements(i),str2double(element_types{i}(end)));
    E_struct.E_ind=zeros(no_elements(i),1);
    E_struct.E_mat=zeros(no_elements(i),1);
    for j=0:1:no_elements(i)-1
        E_struct.E(j+1,:)=sscanf(element_FEB_XML.item(j).getFirstChild.getData.toCharArray()',text_scan_formats{i});
        E_struct.E_ind(j+1)=str2double(element_FEB_XML.item(j).getAttribute('id').toCharArray()');
        E_struct.E_mat(j+1)=str2double(element_FEB_XML.item(j).getAttribute('mat').toCharArray()');
    end    
    E_cell{el_type}=E_struct;
    el_type=el_type+1;
end

%% Alternative method not XML based

% %% Getting element and nodal numbers/coordinates
% 
% %% IMPORT FEB FILE INTO CELL ARRAY
% 
% [T]=txtfile2cell(filename);
% 
% %% Loading .feb as XML
% FEB_XML = xmlread(filename);
% 
% %% Finding material
% MAT_FEB_XML = FEB_XML.getElementsByTagName('Material');
% mat_FEB_XML = MAT_FEB_XML.item(0).getElementsByTagName('material');
% no_mats=mat_FEB_XML.getLength;
% GEO_FEB_XML = FEB_XML.getElementsByTagName('Geometry');
%
% % Finding target fields
% targets={'Material>','Nodes>','Elements>'};
% found_count=[2 2 2];
% [IND_found]=scancell(T,targets,found_count);
% NODE_field=T(IND_found{2}(1):IND_found{2}(2),:);
% ELEMENT_field=T(IND_found{3}(1):IND_found{3}(2),:);
% 
% %Element finding text formats (expand for others)
% text_scan_formats={'<hex8 id="%f" mat="%f">     %f,    %f,    %f,     %f,     %f,    %f,    %f,     %f</hex8>';...
%     '<quad4 id="%f" mat="%f">   %f,   %f,   %f,   %f</quad4>';...
%     '<tet4 id="%f" mat="%f">   %f,   %f,   %f,   %f</tet4>'};
% E={};
% E_num={};
% E_type={};
% for i=1:1:no_mats
%     target=['mat="',num2str(i),'">'];
%     L=cellfun(@(x)~isempty(x),(strfind(ELEMENT_field, target)));
%     if any(L) %if the material has elements
%         element_field=ELEMENT_field(L);
%         element_type=textscan(element_field{1},'%s','delimiter','/>');
%         element_type=element_type{1}{end};
%         
%         switch  element_type
%             case 'hex8'
%                 E_eltype='hex8';
%                 tscan_form=text_scan_formats{1};
%             case 'quad4'
%                 E_eltype='quad4';
%                 tscan_form=text_scan_formats{2};
%             case 'tet4'
%                 E_eltype='tet4';
%                 tscan_form=text_scan_formats{3};
%         end
%         SUB_field=cell2mat(cellfun(@(x) cell2mat(textscan(x,tscan_form)),element_field,'UniformOutput',0));
%         E_elnum=SUB_field(:,1);
%         E_matnum=SUB_field(1,2);
%         E_nodenum=SUB_field(:,3:end);
%         E{i}=E_nodenum;
%         E_num{i}=E_elnum;
%         E_type{i}=E_eltype;
%     else
%         E{i}=[];
%         E_num{i}=[];
%     end
% end
% 
% node_text_scan_format='<node id="%f"> %f, %f, %f</node>';
% N_SUB_field=cell2mat(cellfun(@(x) cell2mat(textscan(x,node_text_scan_format)),NODE_field(2:end-1),'UniformOutput',0));
% N_num=N_SUB_field(:,1);
% N=N_SUB_field(:,2:end);

%%

 
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
