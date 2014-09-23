function [docNode]=addConstraintsLevel_FEB(docNode,FEB_struct)

disp('Adding Constraints level')

rootNode = docNode.getDocumentElement;

constraintNode = docNode.createElement('Constraints');
constraintNode = rootNode.appendChild(constraintNode);

%%

for qC=1:1:numel(FEB_struct.Constraints)
    
    %Set rigid body material ID for current constraint set
    rigidId=FEB_struct.Constraints{qC}.RigidId;
    rigid_node = docNode.createElement('rigid_body'); %create entry
    rigid_node = constraintNode.appendChild(rigid_node); %add entry
    attr = docNode.createAttribute('mat'); %Create attribute
    attr.setNodeValue(sprintf('%u',rigidId)); %Set text
    rigid_node.setAttributeNode(attr); %Add attribute
    
    %Set fix constraints
    if isfield(FEB_struct.Constraints{qC},'Fix');
        for q=1:1:numel(FEB_struct.Constraints{qC}.Fix)
            
            currentBC=FEB_struct.Constraints{qC}.Fix{q}.bc;
            
            fix_node = docNode.createElement('fixed'); %create entry
            fix_node = rigid_node.appendChild(fix_node); %add entry
            
            attr = docNode.createAttribute('bc'); %Create attribute
            attr.setNodeValue(currentBC); %Set text
            fix_node.setAttributeNode(attr); %Add attribute
        end
    end
    
    %Set prescribed constraints
    if isfield(FEB_struct.Constraints{qC},'Prescribe');
        for q=1:1:numel(FEB_struct.Constraints{qC}.Fix)
            
            currentBC=FEB_struct.Constraints{qC}.Prescribe{q}.bc;
            currentLC=FEB_struct.Constraints{qC}.Prescribe{q}.lc;
            currentValue=FEB_struct.Constraints{qC}.Prescribe{q}.Scale;                
                
            prescribed_node = docNode.createElement('prescribed'); %create entry
            prescribed_node = rigid_node.appendChild(prescribed_node); %add entry
            
            attr = docNode.createAttribute('bc'); %Create attribute
            attr.setNodeValue(currentBC); %Set text
            prescribed_node.setAttributeNode(attr); %Add attribute
            
            attr = docNode.createAttribute('lc'); %Create attribute
            attr.setNodeValue(sprintf('%u',currentLC)); %Set text
            prescribed_node.setAttributeNode(attr); %Add attribute
            
            prescribed_node.appendChild(docNode.createTextNode(sprintf('%6.7e',currentValue))); %append data text child
            
        end
    end
    
    %Set force constraints
    if isfield(FEB_struct.Constraints{qC},'Force');
        for q=1:1:numel(FEB_struct.Constraints{qC}.Force)
                        
            currentBC=FEB_struct.Constraints{qC}.Force{q}.bc;

            force_node = docNode.createElement('force'); %create entry
            force_node = rigid_node.appendChild(force_node); %add entry
            
            attr = docNode.createAttribute('bc'); %Create attribute
            attr.setNodeValue(currentBC); %Set text
            force_node.setAttributeNode(attr); %Add attribute
        end
    end
    
  
%     rigidProps=FEB_struct.Constraints{qC}.Properties;
%     rigidPropVals=FEB_struct.Constraints{qC}.Values;
%     rigidType=FEB_struct.Constraints{qC}.Type;
%     loadCurves=FEB_struct.Constraints{qC}.LoadCurveIds;
%     
%     for q=1:1:numel(rigidProps)
%         prop_prop_node = docNode.createElement(rigidProps{q}); %create entry
%         prop_prop_node = prop_node.appendChild(prop_prop_node); %add entry
%         
%         attr = docNode.createAttribute('type'); %Create attribute
%         attr.setNodeValue(rigidType); %Set text
%         prop_prop_node.setAttributeNode(attr); %Add attribute
%         
%         attr = docNode.createAttribute('lc'); %Create attribute
%         attr.setNodeValue(sprintf('%u',loadCurves(q))); %Set text
%         prop_prop_node.setAttributeNode(attr); %Add attribute
%         
%         t_form=repmat('%f, ',1,size(rigidPropVals{q},2)); t_form=t_form(1:end-2);
%         prop_prop_node.appendChild(docNode.createTextNode(sprintf(t_form,rigidPropVals{q}))); %append data text child
%     end
%     
end

end