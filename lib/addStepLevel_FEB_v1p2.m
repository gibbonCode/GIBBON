function [docNode]=addStepLevel_FEB_v1p2(docNode,FEB_struct)

if  isfield(FEB_struct,'Step')
    disp('Adding Step level');
    
    febio_spec = docNode.getDocumentElement;
    
    for q_step=1:1:numel(FEB_struct.Step)
        
        StepNode = docNode.createElement('Step');
        febio_spec.appendChild(StepNode);
                
        attr = docNode.createAttribute('name'); %Create id attribute
        attr.setNodeValue(['Step_',num2str(q_step)]); %Set text
        StepNode.setAttributeNode(attr); %Add attribute
        docRootNode=StepNode;

        %Add module type
        if ~isfield(FEB_struct.Step(q_step),'ModuleType')
            FEB_struct.Step(q_step).ModuleType='solid'; %Use default
        elseif isempty(FEB_struct.Step(q_step).ModuleType)
            FEB_struct.Step(q_step).ModuleType='solid'; %Use default
        end        
        ModuleNode = docNode.createElement('Module'); %create entry
        ModuleNode = docRootNode.appendChild(ModuleNode); %add entry
        attr = docNode.createAttribute('type'); %Create attribute
        attr.setNodeValue(FEB_struct.Step(q_step).ModuleType); %Set text
        ModuleNode.setAttributeNode(attr); %Add id attribute
        
        % Adding boundary components
        if isfield(FEB_struct.Step(q_step),'Boundary')            
            BoundaryNode = docNode.createElement('Boundary'); %create entry
            BoundaryNode = docRootNode.appendChild(BoundaryNode); %add entry
            
            disp('----> Adding Boundary field');
            FEB_struct_temp=FEB_struct.Step(q_step);
            FEB_struct_temp.disp_opt=FEB_struct.disp_opt;
            [docNode]=addBoundaryComponents_FEB(docNode,BoundaryNode,FEB_struct_temp);
        end
      
        % Adding control components
        if isfield(FEB_struct.Step(q_step),'Control')            
            ControlNode = docNode.createElement('Control'); %create entry
            ControlNode = docRootNode.appendChild(ControlNode); %add entry
            
            disp('----> Adding Control field')
            FEB_struct_temp=FEB_struct.Step(q_step);
            FEB_struct_temp.disp_opt=FEB_struct.disp_opt;
            [docNode]=addControlComponents_FEB(docNode,ControlNode,FEB_struct_temp);
        end
        
    end        
end
