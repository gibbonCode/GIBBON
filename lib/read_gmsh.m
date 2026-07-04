

function MeshStruct = read_gmsh(fileName)

T=txtfile2cell(fileName);

start_physicals =  find(contains(T,'$PhysicalNames'));
start_entities = find(contains(T,'$Entities'));
start_nodes = find(contains(T,'$Nodes'));
start_elements = find(contains(T,'$Elements'));

entity_info = str2num(T{start_entities+1});
node_info = str2num(T{start_nodes+1});
element_info = str2num(T{start_elements+1});

num_nodes = node_info(2);
num_node_blocks = node_info(1);
num_elements = element_info(2);
num_element_blocks = element_info(1);
num_entities = sum(entity_info);

%%%read physicals
num_physicals = str2num(T{start_physicals+1});
Physicals = T(start_physicals+2:start_physicals+1+num_physicals);
for i = 1:size(Physicals, 1)
    Physicals{i} = strsplit(Physicals{i},' ');
    Physicals{i}{1} = str2num(Physicals{i}{1});
    Physicals{i}{2} = str2num(Physicals{i}{2});
end

%%%read and sort entities
Entities = cell(4,1);
Points = {};Curves = {};Surfaces = {};Volumes = {};
line = start_entities+2;
for i = 1:entity_info(1)
    Temp = str2num(T{line});
    Points{i,1} = Temp;
    line = line+1;
end
for i = 1:entity_info(2)
    Temp = str2num(T{line});
    Curves{i,1} = Temp;
    line = line+1;
end
for i = 1:entity_info(3)
    Temp =  str2num(T{line});
    Surfaces{i,1} = Temp;
    line = line+1;
end
for i = 1:entity_info(4)
    Temp = str2num(T{line});
    Volumes{i,1} = Temp;
    line = line+1;
end
Entities{1} = Points;
Entities{2} = Curves;
Entities{3} = Surfaces;
Entities{4} = Volumes;

%%%read and sort nodes
node_block_info = zeros(num_node_blocks,4);
node_block_numbers = cell(num_node_blocks,1);
node_block_coords = cell(num_node_blocks,1);
line = start_nodes+2;
for i = 1:num_node_blocks
    
    node_block_info_i = str2num(T{line});
    node_block_size= node_block_info_i(1,4);
    
    node_block_info(i,:) = node_block_info_i;
    node_block_numbers_temp = T(line+1:line+node_block_size);
    node_block_coords_temp = T(line+node_block_size +1:line+node_block_size+node_block_size);
    line = line+2*node_block_size+1;
    
    for j = 1:node_block_size        
        node_block_numbers{i}(j,:) = str2num(node_block_numbers_temp{j});
        node_block_coords{i}(j,:) = str2num(node_block_coords_temp{j});
    end    
       
    
end


%%%Read and sort elements
element_block_info = zeros(num_element_blocks,4);
element_block_numbers = cell(num_element_blocks,1);
line = start_elements+2;
for i = 1:num_element_blocks
    
    element_block_info_i = str2num(T{line});element_block_size= element_block_info_i(1,4);
    element_block_info(i,:) = element_block_info_i;
    element_block_numbers_temp = T(line+1:line+element_block_size);
    
    
    for j = 1:element_block_size        
        element_block_numbers{i}(j,:) = str2num(element_block_numbers_temp{j});
    end    
       
    line = line+element_block_size+1;
end

node_ids = [];
nodes = [];
for i = 1:num_node_blocks  
       
    nodes = [nodes;node_block_coords{i}];
    node_ids = [node_ids;node_block_numbers{i}];
    
end

nodes = double(nodes);
Mesh.Physicals = Physicals;
Mesh.entity_info = entity_info;
Mesh.node_info = node_info;
Mesh.element_info = element_info;
Mesh.nodes = nodes;
Mesh.Node_Blocks = node_block_numbers;
Mesh.Node_ids = node_ids;
Mesh.node_block_info = node_block_info;
Mesh.Entities = Entities;
Mesh.Elements = element_block_numbers;
Mesh.element_block_info = element_block_info; 


%%
vol_ct = 0;
surf_ct = 0;
curve_ct= 0;

[~,part_name,~] = fileparts(fileName);

 node_ids = Mesh.Node_ids; %The node id's
 nodes = Mesh.nodes; %The nodel coordinates

 %%
 
 
for j = 1:size(Mesh.Physicals,1)

    physical = Mesh.Physicals{j};
    phys_name = physical{3}(2:end-1);
    dim = physical{1};
    tag = physical{2};

    switch dim
        case 3

            vol_ct = vol_ct+1;
            entities = Mesh.Entities{4,1};
            element_ids = [];
            element_mat = [];

            for k = 1:size(entities, 1)
                num_phys = entities{k}(8);
                phys_tags = entities{k}(9:9+num_phys-1);
                if sum(find(phys_tags==tag))
                   entity_tag = entities{k}(1);
                   elements_3 = (Mesh.element_block_info(:,1) == dim).*Mesh.element_block_info;
                   entity_ind = find(elements_3(:,2)==entity_tag);
                   element_ids = [element_ids;Mesh.Elements{entity_ind}(:,1)];
                   element_mat = [element_mat;Mesh.Elements{entity_ind}(:,2:end)];                         
                end

            end

            %febio_spec.Mesh.Elements{vol_ct}.ATTR.type='tet4'; %Element type
            volumes_names{vol_ct,1} = phys_name; %Name of this part
            elements_IDs_blocks{vol_ct} = element_ids; %Element id's
            elements_blocks{vol_ct} = element_mat; %The element matrix


        case 2

            surf_ct = surf_ct+1;
            entities = Mesh.Entities{3,1};
            element_mat = []; 
            element_start = 0;

            for k = 1:size(entities, 1)
                num_phys = entities{k}(8);
                phys_tags = entities{k}(9:9+num_phys-1);
                if sum(find(phys_tags==tag))
                   entity_tag = entities{k}(1);
                   elements_2 = (Mesh.element_block_info(:,1) == dim).*(Mesh.element_block_info(:,2) == entity_tag);
                   entity_ind = find(elements_2);
                   element_mat = [element_mat;Mesh.Elements{entity_ind}(:,2:end)];
                   element_start = element_start+1;
                   element_start = element_start+size(element_mat,1);
                end

            end

            element_ids = [1:size(element_mat,1)]';
            facesBoundary_Names{surf_ct,1} = phys_name; %Name of this surface
            facesBoundary_IDs_Block{surf_ct} = element_ids; %Element id's
            facesBoundary_Block{surf_ct} =element_mat; %The element matrix   

            case 1

            curve_ct = curve_ct+1;
            entities = Mesh.Entities{2,1};
            element_mat = []; 
            element_start = 0;

            for k = 1:size(entities, 1)
                num_phys = entities{k}(8);
                phys_tags = entities{k}(9:9+num_phys-1);
                if sum(find(phys_tags==tag))
                   entity_tag = entities{k}(1);
                   elements_1 = (Mesh.element_block_info(:,1) == dim).*(Mesh.element_block_info(:,2) == entity_tag);
                   entity_ind = find(elements_1);
                   element_mat = [element_mat;Mesh.Elements{entity_ind}(:,2:end)];
                   element_start = element_start+1;
                   element_start = element_start+size(element_mat,1);
                end

            end

             element_ids_block{curve_ct} = element_mat(:,1);
             curve_names{curve_ct} = phys_name; %Name of this curve
             curve_elements_block{curve_ct} = element_ids; %Element id's

   end
end

%%
facesBoundary = [];
boundaryMarker = [];
elements = [];
elementsMaterialID = [];

%%
for i = 1:vol_ct
    
    elementMaterialID_block = ones(length(elements_IDs_blocks{i}),1).*i;
    elements = [elements;elements_blocks{i}];
    elementsMaterialID = [elementsMaterialID;elementMaterialID_block];
    
end

for i = 1:surf_ct
    
    boundaryMarker_block = ones(length(facesBoundary_IDs_Block{i}),1).*i;
    facesBoundary = [facesBoundary;facesBoundary_Block{i}];
    boundaryMarker = [boundaryMarker;boundaryMarker_block];
    
end

% for i = 1:curve_ct
%     
% end

MeshStruct.node_IDs = node_ids;
MeshStruct.nodes = nodes;
MeshStruct.facesBoundary = facesBoundary;
MeshStruct.boundaryMarker = boundaryMarker;
MeshStruct.elements = elements;
MeshStruct.elementMaterialID = elementsMaterialID;
MeshStruct.loadNameStruct.MeshName = part_name;
MeshStruct.loadNameStruct.VolumeNames = volumes_names;
MeshStruct.loadNameStruct.SurfaceNames = facesBoundary_Names;


%                 nodes: [415×3 double]
%         facesBoundary: [320×3 double]
%        boundaryMarker: [320×1 double]
%                 faces: [8144×3 double]
%              elements: [2036×4 double]
%     elementMaterialID: [2036×1 double]
%        faceMaterialID: [8144×1 double]
%        loadNameStruct: [1×1 struct]

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
