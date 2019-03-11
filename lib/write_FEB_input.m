function docNode=write_FEB_input(FEB_struct)

% function docNode=write_FEB_input(FEB_struct)
% ------------------------------------------------------------------------
%
% This function uses the input FEB_struct to generate a basic .feb
% file suitable for import into PreView. 
%
% %% EXAMPLE
% This example generates a 4 element model (2 hex8, 1 quad4and 1 tet4),
% where each element has a different material.
%
% savename='test.feb';%
% N=[0.5,0.5,1.5; 0.5,1.5,1.5; 1.5,1.5,1.5; 1.5,0.5,1.5; 0.5,0.5,0.5; 0.5,1.5,0.5; 1.5,1.5,0.5; 1.5,0.5,0.5; 0.5,2.5,1.5; 1.5,2.5,1.5; 0.5,2.5,0.5; 1.5,2.5,0.5;];
% E={[1 2 3 4 5 6 7 8; 2 9 10 3 6 11 12 7;]; [1 2 3 4]; [1 2 3 5];};
% E_type={'hex8';'quad4';'tet4';};
% M={[1 2];3;4};
% FEB_struct.Geometry.Nodes=N;
% FEB_struct.Geometry.Elements=E;
% FEB_struct.Geometry.ElementType=E_type;
% FEB_struct.Geometry.ElementMat=M;
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 10/08/2010
%------------------------------------------------------------------------
% 
%%

%% Initialize docNode object
docNode = com.mathworks.xml.XMLUtils.createDocument('febio_spec'); %Create the overall febio_spec field
febio_spec = docNode.getDocumentElement;
febio_spec.setAttribute('version','1.2'); %Addinng version attribut (important for interpretation)

%% DEFINING MATERIAL LEVEL
ElementAddNode = docNode.createElement('Material');
febio_spec.appendChild(ElementAddNode);
MATNode = docNode.getElementsByTagName('Material').item(0);

% Adding material fields
disp('Adding Material level')
A=[FEB_struct.Geometry.ElementMat{:}];
mat_ind=unique(A);
parent_node =MATNode;

hw = waitbar(0,'Adding material level...'); 

for i=1:1:numel(mat_ind)
    
    %Adding material node
    material_node = docNode.createElement('material'); %create material entry
    material_node = parent_node.appendChild(material_node); %add material entry
    
    %Checking material type
    Mat_i=FEB_struct.Materials{i};
    mat_type=Mat_i.type;
    switch Mat_i.type
        case {'uncoupled solid mixture','solid mixture'}
            
            %Mat ID
            attr = docNode.createAttribute('id'); %Create attribute
            attr.setNodeValue(num2str(mat_ind(i))); %Set text
            material_node.setAttributeNode(attr); %Add attribute
            
            %Mat Name
            attr = docNode.createAttribute('name'); %Create attribute
            attr.setNodeValue(['mat_',num2str(i)]); %Set text
            material_node.setAttributeNode(attr); %Add attribute
            
            %Mat Type
            attr = docNode.createAttribute('type'); %Create attribute
            attr.setNodeValue(mat_type); %Set text
            material_node.setAttributeNode(attr); %Add attribute
            
            for j=1:1:numel(FEB_struct.Materials{i}.Mats)
                
                solid_node = docNode.createElement('solid'); %create entry
                solid_node = material_node.appendChild(solid_node); %add entry
                
                %Mat Type
                attr = docNode.createAttribute('type'); %Create attribute
                attr.setNodeValue(Mat_i.Mats{j}.type); %Set text
                solid_node.setAttributeNode(attr); %Add attribute
                
                %Mat ID
                attr = docNode.createAttribute('id'); %Create attribute
                attr.setNodeValue(num2str(j)); %Set text
                solid_node.setAttributeNode(attr); %Add attribute
                
                mat_props=Mat_i.Mats{j}.props;
                mat_prop_vals=Mat_i.Mats{j}.vals;
                %Material properties
                for p=1:1:numel(mat_props)
                    mat_prop=mat_props{p};
                    mat_prop_val=mat_prop_vals{p};
                    prop_node = docNode.createElement(mat_prop); %create entry
                    prop_node = solid_node.appendChild(prop_node); %add entry
                    t_form=repmat('%6.7e, ',1,size(mat_prop_val,2)); t_form=t_form(1:end-2);
                    prop_node.appendChild(docNode.createTextNode(sprintf(t_form,mat_prop_val))); %append data text child
                end  
                
                %Fiber spec
                mat_aniso_type=Mat_i.Mats{j}.aniso_type;
                switch mat_aniso_type
                    case 'none'
                        
                    case 'fiber'
                        prop_node = docNode.createElement('fiber'); %create entry
                        prop_node = solid_node.appendChild(prop_node); %add entry
                        attr = docNode.createAttribute('type'); %Create attribute
                        attr.setNodeValue('user'); %Set text
                        prop_node.setAttributeNode(attr); %Add attribute
                    case 'mat_axis'
                        prop_node = docNode.createElement('mat_axis'); %create entry
                        prop_node = solid_node.appendChild(prop_node); %add entry
                        attr = docNode.createAttribute('type'); %Create attribute
                        attr.setNodeValue('vector'); %Set text
                        prop_node.setAttributeNode(attr); %Add attribute
                end
                
            end
            
        otherwise            
            
            mat_props=Mat_i.props;
            mat_prop_vals=Mat_i.vals;
            mat_aniso_type=Mat_i.aniso_type;            
            
            %Mat ID
            attr = docNode.createAttribute('id'); %Create attribute
            attr.setNodeValue(num2str(mat_ind(i))); %Set text
            material_node.setAttributeNode(attr); %Add attribute
            
            %Mat Name
            attr = docNode.createAttribute('name'); %Create attribute
            attr.setNodeValue(['mat_',num2str(i)]); %Set text
            material_node.setAttributeNode(attr); %Add attribute
            
            %Mat Type
            attr = docNode.createAttribute('type'); %Create attribute
            attr.setNodeValue(mat_type); %Set text
            material_node.setAttributeNode(attr); %Add attribute
            
            %Material properties
            for p=1:1:numel(mat_props)
                mat_prop=mat_props{p};
                mat_prop_val=mat_prop_vals{p};
                prop_node = docNode.createElement(mat_prop); %create entry
                prop_node = material_node.appendChild(prop_node); %add entry
                t_form=repmat('%6.7e, ',1,size(mat_prop_val,2)); t_form=t_form(1:end-2);
                prop_node.appendChild(docNode.createTextNode(sprintf(t_form,mat_prop_val))); %append data text child
            end
            
            %Fiber spec
            switch mat_aniso_type
                case 'none' %No fibre type specified, do nothing
                    
                case 'fiber'
                    prop_node = docNode.createElement('fiber'); %create entry
                    prop_node = material_node.appendChild(prop_node); %add entry
                    attr = docNode.createAttribute('type'); %Create attribute
                    attr.setNodeValue('user'); %Set text
                    prop_node.setAttributeNode(attr); %Add attribute
                case 'mat_axis'
                    prop_node = docNode.createElement('mat_axis'); %create entry
                    prop_node = material_node.appendChild(prop_node); %add entry
                    attr = docNode.createAttribute('type'); %Create attribute
                    attr.setNodeValue('vector'); %Set text
                    prop_node.setAttributeNode(attr); %Add attribute
            end
    end
    %     prop_node.appendChild(docNode.createTextNode(' ')); %append data text child
    waitbar(i/numel(mat_ind),hw,'Adding material level...'); 
end
close(hw);
drawnow; 

%% DEFINING GEOMETRY LEVEL
disp('Adding Geometry level')
ElementAddNode = docNode.createElement('Geometry');
febio_spec.appendChild(ElementAddNode);
GEONode = docNode.getElementsByTagName('Geometry').item(0);

% Adding Node field
disp('----> Adding node field')

parent_node = docNode.createElement('Nodes');
parent_node = GEONode.appendChild(parent_node);

hw = waitbar(0,'Adding node entries....');
n_steps=size(FEB_struct.Geometry.Nodes,1);
for q_n=1:1:n_steps
    node_node = docNode.createElement('node'); %create node entry
    node_node = parent_node.appendChild(node_node); %add node entry
    attr = docNode.createAttribute('id'); %Create id attribute
    attr.setNodeValue(num2str(q_n)); %Set id text
    node_node.setAttributeNode(attr); %Add id attribute
    node_node.appendChild(docNode.createTextNode(sprintf('%6.7e, %6.7e, %6.7e',FEB_struct.Geometry.Nodes(q_n,:)))); %append data text child
    waitbar(q_n/n_steps);
end
close(hw);

% Adding Elements field
disp('----> Adding element field')

parent_node = docNode.createElement('Elements');
parent_node = GEONode.appendChild(parent_node);
LogicSetThickness=0;
if isfield(FEB_struct.Geometry,'ElementData');
    if isfield(FEB_struct.Geometry.ElementData,'Thickness');
        LogicSetThickness=1;
    end
end

addedElementDataField=0;
e_ind=1;
for q_e1=1:1:numel(FEB_struct.Geometry.Elements)
    E=FEB_struct.Geometry.Elements{q_e1};
    E_type=FEB_struct.Geometry.ElementType{q_e1};
    M=FEB_struct.Geometry.ElementMat{q_e1};
    
    hw = waitbar(0,['Adding ',E_type,' element entries....']);
    disp(['----> Adding ',E_type,' element entries....'])
    n_steps=size(E,1);
    for q_e2=1:1:n_steps
        element_node = docNode.createElement(E_type); %create element entry
        element_node = parent_node.appendChild(element_node); %add element entry
        
        attr = docNode.createAttribute('id'); %Create id attribute
        attr.setNodeValue(num2str(e_ind)); %Set id text
        element_node.setAttributeNode(attr); %Add id attribute
        
        attr = docNode.createAttribute('mat'); %Create mat attribute
        attr.setNodeValue(num2str(M(q_e2))); %Set mat text
        element_node.setAttributeNode(attr); %Add mat attribute
        
        t_form=repmat('   %i,',1,size(E,2)); t_form=t_form(1:end-1);
        element_node.appendChild(docNode.createTextNode(sprintf(t_form,E(q_e2,:)))); %append data text child
        
        %Adding element data entries for quad4 elements
        if strcmp(E_type,'quad4') && LogicSetThickness            
            if ismember(e_ind,FEB_struct.Geometry.ElementData.IndicesForThickness)
                %Get Geometry level
                docRootNode =  GEONode;
                if addedElementDataField==0
                    parent_node1 = docNode.createElement('ElementData');
                    parent_node1 = docRootNode.appendChild(parent_node1);
                    addedElementDataField=1;
                end
                
                %Creating ElementData entries
                element_node = docNode.createElement('element');
                element_node = parent_node1.appendChild(element_node);
                attr = docNode.createAttribute('id');
                attr.setNodeValue(num2str(e_ind));
                element_node.setAttributeNode(attr);
                
                %Adding thickness level
                v_text=sprintf('%6.7e, %6.7e, %6.7e, %6.7e',FEB_struct.Geometry.ElementData.Thickness(FEB_struct.Geometry.ElementData.IndicesForThickness==e_ind)*ones(1,4));
                thickness_node = docNode.createElement('thickness');
                element_node.appendChild(thickness_node);
                thickness_node.appendChild(docNode.createTextNode(v_text)); %Adding thickness data text
            end
        end
   
        %Adding element data entries for tri3 elements
        if strcmp(E_type,'tri3') && LogicSetThickness            
            if ismember(e_ind,FEB_struct.Geometry.ElementData.IndicesForThickness)
                %Get Geometry level
                docRootNode =  GEONode;
                if addedElementDataField==0
                    parent_node1 = docNode.createElement('ElementData');
                    parent_node1 = docRootNode.appendChild(parent_node1);
                    addedElementDataField=1;
                end
                
                %Creating ElementData entries
                element_node = docNode.createElement('element');
                element_node = parent_node1.appendChild(element_node);
                attr = docNode.createAttribute('id');
                attr.setNodeValue(num2str(e_ind));
                element_node.setAttributeNode(attr);
                
                %Adding thickness level
                v_text=sprintf('%6.7e, %6.7e, %6.7e',FEB_struct.Geometry.ElementData.Thickness(FEB_struct.Geometry.ElementData.IndicesForThickness==e_ind)*ones(1,3));
                thickness_node = docNode.createElement('thickness');
                element_node.appendChild(thickness_node);
                thickness_node.appendChild(docNode.createTextNode(v_text)); %Adding thickness data text
            end
        end
        
        %%
        e_ind=e_ind+1;
        waitbar(q_e2/n_steps);
    end
    close(hw);
end

% Checking for fiber / mat_axis information
if isfield(FEB_struct.Geometry,'ElementData') %If ElementData exists
    if isfield(FEB_struct.Geometry.ElementData,'MatAxis') %If MatAxis data exists
        docNode=add_fiber_dir_FEB(docNode,FEB_struct,0); %Adding fibre directions
    end
end

%% Adding Output level
disp('Adding Output level')

ElementAddNode = docNode.createElement('Output');
febio_spec.appendChild(ElementAddNode);
OUTPUTNode = docNode.getElementsByTagName('Output').item(0);

%% Adding Plotfile field
disp('----> Adding Plotfile field')
parent_node = docNode.createElement('plotfile');
parent_node = OUTPUTNode.appendChild(parent_node);
parent_node.setAttribute('type','febio');

% Adding Output requests
if isfield(FEB_struct,'Output');
    for q=1:1:numel(FEB_struct.Output.VarTypes)
        var_node = docNode.createElement('var'); %create entry
        var_node = parent_node.appendChild(var_node); %add entry
        var_node.setAttribute('type',FEB_struct.Output.VarTypes{q}); %add attribute
    end
end

%% Checking for contact information
disp('Boundary condition level')
if isfield(FEB_struct,'Boundary') %If boundary field exists
    if isfield(FEB_struct.Boundary,'Contact') %If MatAxis data exists
        docNode=add_contact_pair_FEB(docNode,FEB_struct,0); %Adding contact
    end
end


%% 

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
