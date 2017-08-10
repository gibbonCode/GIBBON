function [domNode]=addConstraintsLevel_FEB(domNode,FEB_struct)

% function [domNode]=addConstraintsLevel_FEB(domNode,FEB_struct)
% ------------------------------------------------------------------------
% Adds contraints information to the XML object domNode based on
% the FEBio structure FEB_struct. 
% 
% Fixed and prescribed constraints are supported
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2016/06/02
%------------------------------------------------------------------------

%%
disp('Adding Constraints level')

rootNode = domNode.getDocumentElement;

constraintNode = domNode.createElement('Constraints');
constraintNode = rootNode.appendChild(constraintNode);

%%

for qC=1:1:numel(FEB_struct.Constraints)
    
    %Set rigid body material ID for current constraint set
    rigidId=FEB_struct.Constraints{qC}.RigidId;
    rigid_node = domNode.createElement('rigid_body'); %create entry
    rigid_node = constraintNode.appendChild(rigid_node); %add entry
    attr = domNode.createAttribute('mat'); %Create attribute
    attr.setNodeValue(sprintf('%u',rigidId)); %Set text
    rigid_node.setAttributeNode(attr); %Add attribute
    
    %Set fix constraints
    if isfield(FEB_struct.Constraints{qC},'Fix');
        for q=1:1:numel(FEB_struct.Constraints{qC}.Fix)
            
            currentBC=FEB_struct.Constraints{qC}.Fix{q}.bc;
            
            fix_node = domNode.createElement('fixed'); %create entry
            fix_node = rigid_node.appendChild(fix_node); %add entry
            
            attr = domNode.createAttribute('bc'); %Create attribute
            attr.setNodeValue(currentBC); %Set text
            fix_node.setAttributeNode(attr); %Add attribute
        end
    end
    
    %Set prescribed constraints
    if isfield(FEB_struct.Constraints{qC},'Prescribe');
        for q=1:1:numel(FEB_struct.Constraints{qC}.Prescribe)
            
            currentBC=FEB_struct.Constraints{qC}.Prescribe{q}.bc;
            currentLC=FEB_struct.Constraints{qC}.Prescribe{q}.lc;
            currentValue=FEB_struct.Constraints{qC}.Prescribe{q}.Scale;                
                
            prescribed_node = domNode.createElement('prescribed'); %create entry
            prescribed_node = rigid_node.appendChild(prescribed_node); %add entry
            
            attr = domNode.createAttribute('bc'); %Create attribute
            attr.setNodeValue(currentBC); %Set text
            prescribed_node.setAttributeNode(attr); %Add attribute
            
            attr = domNode.createAttribute('lc'); %Create attribute
            attr.setNodeValue(sprintf('%u',currentLC)); %Set text
            prescribed_node.setAttributeNode(attr); %Add attribute
            
            prescribed_node.appendChild(domNode.createTextNode(sprintf('%6.7e',currentValue))); %append data text child
            
        end
    end
    
    %Set force constraints
    if isfield(FEB_struct.Constraints{qC},'Force');
        for q=1:1:numel(FEB_struct.Constraints{qC}.Force)
                        
            currentBC=FEB_struct.Constraints{qC}.Force{q}.bc;
            currentLC=FEB_struct.Constraints{qC}.Force{q}.lc;
            currentValue=FEB_struct.Constraints{qC}.Force{q}.Scale;

            force_node = domNode.createElement('force'); %create entry
            force_node = rigid_node.appendChild(force_node); %add entry
            
            attr = domNode.createAttribute('bc'); %Create attribute
            attr.setNodeValue(currentBC); %Set text
            force_node.setAttributeNode(attr); %Add attribute
            
            attr = domNode.createAttribute('lc'); %Create attribute
            attr.setNodeValue(sprintf('%u',currentLC)); %Set text
            force_node.setAttributeNode(attr); %Add attribute
            
            force_node.appendChild(domNode.createTextNode(sprintf('%6.7e',currentValue))); %append data text child
        end
    end
    
  
%     rigidProps=FEB_struct.Constraints{qC}.Properties;
%     rigidPropVals=FEB_struct.Constraints{qC}.Values;
%     rigidType=FEB_struct.Constraints{qC}.Type;
%     loadCurves=FEB_struct.Constraints{qC}.LoadCurveIds;
%     
%     for q=1:1:numel(rigidProps)
%         prop_prop_node = domNode.createElement(rigidProps{q}); %create entry
%         prop_prop_node = prop_node.appendChild(prop_prop_node); %add entry
%         
%         attr = domNode.createAttribute('type'); %Create attribute
%         attr.setNodeValue(rigidType); %Set text
%         prop_prop_node.setAttributeNode(attr); %Add attribute
%         
%         attr = domNode.createAttribute('lc'); %Create attribute
%         attr.setNodeValue(sprintf('%u',loadCurves(q))); %Set text
%         prop_prop_node.setAttributeNode(attr); %Add attribute
%         
%         t_form=repmat('%f, ',1,size(rigidPropVals{q},2)); t_form=t_form(1:end-2);
%         prop_prop_node.appendChild(domNode.createTextNode(sprintf(t_form,rigidPropVals{q}))); %append data text child
%     end
%     
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
