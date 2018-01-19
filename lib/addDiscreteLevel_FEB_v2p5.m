function domNode=addDiscreteLevel_FEB(domNode,FEB_struct)

% function [domNode]=addDiscreteLevel_FEB(domNode,FEB_struct)
% ------------------------------------------------------------------------
% Adds geometry information to the XML object domNode based on
% the FEBio structure FEB_struct.
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2016/06/02
%------------------------------------------------------------------------

%% Set default display setting if missing
if ~isfield(FEB_struct,'disp_opt')
    FEB_struct.disp_opt=0;
end

%% Adding discrete_discrete field

% <Discrete>
%
%   <discrete_discrete id="1" type="nonlinear spring">
%       <force lc="1">1.0</force>
%   </discrete_discrete>
%
%   <discrete_discrete id="2" type="linear spring">
%       <E>2.0</E>
%   </discrete_discrete>
%
%   <discrete dmat="1" discrete_set="springs1"/>
%
%   <discrete dmat="2" discrete_set="springs2">
%
% </Discrete>

%% Define discrete materials

rootNode = domNode.getDocumentElement;

DiscreteNode = domNode.createElement('Discrete');
DiscreteNode = rootNode.appendChild(DiscreteNode);

% Adding discrete fields
disp('Adding Discrete level')

if FEB_struct.disp_opt==1
    hw = waitbar(0,'Adding discrete level...');
end

for q_mat=1:1:numel(FEB_struct.Discrete.discrete_material)
    
    %Adding discrete_material node
    discrete_node = domNode.createElement('discrete_material'); %create discrete entry
    discrete_node = DiscreteNode.appendChild(discrete_node); %add discrete entry
    
    %The current discrete_material struct
    currentDiscreteStruct=FEB_struct.Discrete.discrete_material{q_mat};
    
    %Current discrete name
%     if ~isfield(currentDiscreteStruct,'Name')
%         currentDiscreteStruct.Name=['mat_',sprintf('%u',q_mat)];
%     end
    
    %Current id
    currentDiscreteStruct.id=q_mat;
    
    [domNode]=setDiscreteEntries(domNode,discrete_node,currentDiscreteStruct);
       
    if FEB_struct.disp_opt==1
        %Adjust waitbar
        waitbar(q_mat/numel(matIndicesUnique),hw,'Adding discrete level...');
    end
    
end

if FEB_struct.disp_opt==1
    %Close waitbat
    close(hw);
    drawnow;
end

%% Assign material id's for discrete sets

for q_set=1:1:numel(FEB_struct.Geometry.DiscreteSet)
    
    currentDiscreteStruct=FEB_struct.Geometry.DiscreteSet{q_set};
    
    %Adding discrete node
    discrete_node = domNode.createElement('discrete'); %create discrete entry
    discrete_node = DiscreteNode.appendChild(discrete_node); %add discrete entry
    
    %Setting dmat ID
    attr = domNode.createAttribute('dmat'); %Create attribute
    attr.setNodeValue(sprintf('%u',currentDiscreteStruct.dmat)); %Set text
    discrete_node.setAttributeNode(attr); %Add attribute
    
    %Setting discrete_set
    attr = domNode.createAttribute('discrete_set'); %Create attribute
    attr.setNodeValue(currentDiscreteStruct.Name); %Set text
    discrete_node.setAttributeNode(attr); %Add attribute
end

end

function [domNode]=setDiscreteEntries(domNode,levelNode,inputStruct)

%% ATTRIBUTES
if isfield(inputStruct,'id')
    %Setting discrete ID
    attr = domNode.createAttribute('id'); %Create attribute
    attr.setNodeValue(sprintf('%u',inputStruct.id)); %Set text
    levelNode.setAttributeNode(attr); %Add attribute
end

if isfield(inputStruct,'Type')
    %Setting discrete Type
    attr = domNode.createAttribute('type'); %Create attribute
    attr.setNodeValue(inputStruct.Type); %Set text
    levelNode.setAttributeNode(attr); %Add attribute
end

if isfield(inputStruct,'Name')
    %Setting discrete Name
    attr = domNode.createAttribute('name'); %Create attribute
    attr.setNodeValue(inputStruct.Name); %Set text
    levelNode.setAttributeNode(attr); %Add attribute
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
        prop_node = domNode.createElement(currentProp); %create entry
        prop_node = levelNode.appendChild(prop_node); %add entry
        t_form=repmat('%6.7e, ',1,size(currentVal,2)); t_form=t_form(1:end-2);
        prop_node.appendChild(domNode.createTextNode(sprintf(t_form,currentVal))); %append data text child
        
        %Add potential property attributes (e.g. property load curves)
        if ~isempty(mat_prop_attr_name)
            if q<=numel(mat_prop_attr_name)
                currentPropAttrName=mat_prop_attr_name{q};
                if ~isempty(currentPropAttrName)
                    currentPropAttrVal=mat_prop_attr_val{q};
                    
                    attr = domNode.createAttribute(currentPropAttrName); %Create attribute
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
                    prop_node_sub = domNode.createElement(currentProp); %create entry
                    prop_node_sub = prop_node.appendChild(prop_node_sub); %add entry
                    t_form=repmat('%6.7e, ',1,size(currentVal,2)); t_form=t_form(1:end-2);
                    prop_node_sub.appendChild(domNode.createTextNode(sprintf(t_form,currentVal))); %append data text child
                end
            end
        end
        
    end
end

end
 
%% <-- GIBBON footer text -->
%
%     GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
%     image segmentation, image-based modeling, meshing, and finite element
%     analysis.
%
%     Copyright (C) 2017  Kevin Mattheus Moerman
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2017  Kevin Mattheus Moerman
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
