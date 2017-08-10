function docNode=addMaterialLevel_FEB(docNode,FEB_struct)

% function docNode=addMaterialLevel_FEB(docNode,FEB_struct)
% ------------------------------------------------------------------------
%
% This function uses the input FEB_struct to add the material section to a
% docNode for an FEBio .feb file.
%
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% Change log:
% 2015/06/16: Altered the way viscoelasticity is defined to accomodate
% viscoelasticity for solid mixtures
%------------------------------------------------------------------------

%%

rootNode = docNode.getDocumentElement;

MaterialNode = docNode.createElement('Material');
MaterialNode = rootNode.appendChild(MaterialNode);

% Adding material fields
disp('Adding Material level')
materialIndicesUniquePerSet=cellfun(@(x) unique(x(:))',FEB_struct.Geometry.ElementMat,'UniformOutput',0);
matIndicesUnique=unique([materialIndicesUniquePerSet{:}]);
% materialIndexSet=[FEB_struct.Geometry.ElementMat{:}];
% matIndicesUnique=unique(materialIndexSet);

if FEB_struct.disp_opt==1
    hw = waitbar(0,'Adding material level...');
end

for q_mat=1:1:numel(matIndicesUnique)
    
    %Adding material node
    material_node = docNode.createElement('material'); %create material entry
    material_node = MaterialNode.appendChild(material_node); %add material entry
    
    %The current material struct
    currentMaterialStruct=FEB_struct.Materials{q_mat};
    
    %Current material name
    if ~isfield(currentMaterialStruct,'Name')
        currentMaterialStruct.Name=['mat_',sprintf('%u',q_mat)];
    end
    
    %Current id
    currentMaterialStruct.id=matIndicesUnique(q_mat);
    
    [docNode]=setMaterialEntries(docNode,material_node,currentMaterialStruct);
    
    %% Adding fields for solid mixture or multigeneration materials
    
    switch currentMaterialStruct.Type
        case {'uncoupled solid mixture','solid mixture','biphasic'}
            numSolids=numel(currentMaterialStruct.Solid);
            for q_sol=1:1:numSolids
                %Add solid level
                solid_node = docNode.createElement('solid'); %create entry
                solid_node = material_node.appendChild(solid_node); %add entry
                
                %Set parameters
                currentSolidStruct=currentMaterialStruct.Solid{q_sol};
                [docNode]=setMaterialEntries(docNode,solid_node,currentSolidStruct);
            end
        case 'multigeneration'
            numGenerations=numel(currentMaterialStruct.Generation);
            for q_gen=1:1:numGenerations
                currentGenerationStruct=currentMaterialStruct.Generation{q_gen};
                
                %Current id
                currentGenerationStruct.id=q_gen;
                
                %Add solid level
                gen_node = docNode.createElement('generation'); %create entry
                gen_node = material_node.appendChild(gen_node); %add entry
                
                [docNode]=setMaterialEntries(docNode,gen_node,currentGenerationStruct);
                
                numSolids=numel(currentMaterialStruct.Generation{q_gen}.Solid);
                for q_sol=1:1:numSolids
                    %Add solid level
                    solid_node = docNode.createElement('solid'); %create entry
                    solid_node = gen_node.appendChild(solid_node); %add entry
                    
                    %Set parameters
                    currentSolidStruct=currentGenerationStruct.Solid{q_sol};
                    [docNode]=setMaterialEntries(docNode,solid_node,currentSolidStruct);
                end
            end
        case {'uncoupled viscoelastic','viscoelastic'}
            if isfield(currentMaterialStruct,'Elastic')
                %Add elastic level
                elastic_node = docNode.createElement('elastic'); %create entry
                elastic_node = material_node.appendChild(elastic_node); %add entry
                
                %Get material structure for elastic
                elasticStruct=currentMaterialStruct.Elastic;
                
                %Setting elastic material Type
                attr = docNode.createAttribute('type'); %Create attribute
                attr.setNodeValue(elasticStruct.Type); %Set text
                elastic_node.setAttributeNode(attr); %Add attribute
                
                switch elasticStruct.Type
                    case {'uncoupled solid mixture','solid mixture'}
                        numSolids=numel(elasticStruct.Solid);
                        for q_sol=1:1:numSolids
                            %Add solid level
                            solid_node = docNode.createElement('solid'); %create entry
                            solid_node = elastic_node.appendChild(solid_node); %add entry
                            
                            %Set parameters
                            currentSolidStruct=elasticStruct.Solid{q_sol};
                            [docNode]=setMaterialEntries(docNode,solid_node,currentSolidStruct);
                        end
                    otherwise
                        [docNode]=setMaterialEntries(docNode,elastic_node,elasticStruct);
                end
            else
                warning('Elastic field (containing material structure) missing. Assuming "old" viscoelastic specification type. It is recommended to update code.')
            end
    end
    
    if FEB_struct.disp_opt==1
        %Adjust waitbar
        waitbar(q_mat/numel(matIndicesUnique),hw,'Adding material level...');
    end
    
end

if FEB_struct.disp_opt==1
    %Close waitbat
    close(hw);
    drawnow;
end

end


function [docNode]=setMaterialEntries(docNode,levelNode,inputStruct)

%% ATTRIBUTES
if isfield(inputStruct,'id')
    %Setting material ID
    attr = docNode.createAttribute('id'); %Create attribute
    attr.setNodeValue(sprintf('%u',inputStruct.id)); %Set text
    levelNode.setAttributeNode(attr); %Add attribute
end

if isfield(inputStruct,'Type')
    %Setting material Type
    attr = docNode.createAttribute('type'); %Create attribute
    attr.setNodeValue(inputStruct.Type); %Set text
    levelNode.setAttributeNode(attr); %Add attribute
end

if isfield(inputStruct,'Name')
    %Setting material Name
    attr = docNode.createAttribute('name'); %Create attribute
    attr.setNodeValue(inputStruct.Name); %Set text
    levelNode.setAttributeNode(attr); %Add attribute
end

%% ANISOTROPY PARAMETERS

if isfield(inputStruct,'AnisoType')
    mat_aniso_type=inputStruct.AnisoType;
else %DEFAULT
    mat_aniso_type='none';
end

%Fiber spec
switch mat_aniso_type
    case 'none' %No fibre type specified, do nothing
        
    case 'fiber'
        prop_node = docNode.createElement('fiber'); %create entry
        prop_node = levelNode.appendChild(prop_node); %add entry
        attr = docNode.createAttribute('type'); %Create attribute
        attr.setNodeValue('user'); %Set text
        prop_node.setAttributeNode(attr); %Add attribute
    case 'mat_axis'
        prop_node = docNode.createElement('mat_axis'); %create entry
        prop_node = levelNode.appendChild(prop_node); %add entry
        attr = docNode.createAttribute('type'); %Create attribute
        attr.setNodeValue('user'); %Set text
        prop_node.setAttributeNode(attr); %Add attribute
end

%% PROPERTIES
if isfield(inputStruct,'Properties')
    
    %Access properties and values
    mat_props=inputStruct.Properties;
    mat_prop_vals=inputStruct.Values;
    
    %Access property attributes
    if isfield(inputStruct,'PropAttrName')
        mat_prop_attr_name=inputStruct.PropAttrName;
        mat_prop_attr_val=inputStruct.PropAttrVal;
    else
        mat_prop_attr_name=[];
        mat_prop_attr_val=[];
    end
    
    %Access property parameters
    if isfield(inputStruct,'PropParName')
        mat_prop_par_name=inputStruct.PropParName;
        mat_prop_par_val=inputStruct.PropParVal;
    else
        mat_prop_par_name=cell(1,numel(mat_props));
        mat_prop_par_val=cell(1,numel(mat_props));
    end
    
    %Set property values
    for q=1:1:numel(mat_props)
        currentProp=mat_props{q};
        currentVal=mat_prop_vals{q};
        prop_node = docNode.createElement(currentProp); %create entry
        prop_node = levelNode.appendChild(prop_node); %add entry
        t_form=repmat('%6.7e, ',1,size(currentVal,2)); t_form=t_form(1:end-2);
        prop_node.appendChild(docNode.createTextNode(sprintf(t_form,currentVal))); %append data text child
        
        %Add potential property attributes (e.g. property load curves)
        if ~isempty(mat_prop_attr_name)
            if q<=numel(mat_prop_attr_name)
                currentPropAttrName=mat_prop_attr_name{q};
                if ~isempty(currentPropAttrName)
                    currentPropAttrVal=mat_prop_attr_val{q};
                    
                    attr = docNode.createAttribute(currentPropAttrName); %Create attribute
                    if isfloat(currentPropAttrVal)
                        t_form=repmat('%6.7e, ',1,size(currentPropAttrVal,2));
                        t_form=t_form(1:end-2);
                        strEntry=sprintf(t_form,currentPropAttrVal);
                    elseif isa(currentPropAttrVal,'char')
                        strEntry=currentPropAttrVal;
                    else
                        error('Unknown class for currentPropAttrVal');
                    end
                    attr.setNodeValue(strEntry); %Set text
                    prop_node.setAttributeNode(attr); %Add attribute
                end
            end
        end
        
        %Add potential property parameters
        if ~isempty(mat_prop_par_name{q})
            for qp=1:1:numel(mat_prop_par_name{q})
                currentProp=mat_prop_par_name{q}{qp};
                currentVal=mat_prop_par_val{q}{qp};
                if ~isempty(currentProp)
                    prop_node_sub = docNode.createElement(currentProp); %create entry
                    prop_node_sub = prop_node.appendChild(prop_node_sub); %add entry
                    t_form=repmat('%6.7e, ',1,size(currentVal,2)); t_form=t_form(1:end-2);
                    prop_node_sub.appendChild(docNode.createTextNode(sprintf(t_form,currentVal))); %append data text child
                end
            end
        end
        
    end
end

end
 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
