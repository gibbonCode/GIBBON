function [docNode]=addBoundaryComponents_FEB_v2p5(docNode,boundaryNode,FEB_struct)


%% ADDING FIX BC's

if  isfield(FEB_struct.Boundary,'Fix')
    disp('----> Defining fix type boundary conditions');
    
    for q=1:1:numel(FEB_struct.Boundary.Fix) %For all the fix type BC's
        
        currentBC=FEB_struct.Boundary.Fix{q}.bc;
        
        %Create fix section
        parent_node = docNode.createElement('fix');
        parent_node = boundaryNode.appendChild(parent_node);
        
        %Set bc attribute
        attr = docNode.createAttribute('bc'); %Create id attribute
        attr.setNodeValue(currentBC); %Set id text
        parent_node.setAttributeNode(attr); %Add id attribute
        
        %Check consistency
        if isfield(FEB_struct.Boundary.Fix{q},'SetName') && isfield(FEB_struct.Boundary.Fix{q},'Set')
            error('Specify either SetName or Set, not both');
        end
        
        %Add set name attribute if present
        if isfield(FEB_struct.Boundary.Fix{q},'SetName')
            currentSetName=FEB_struct.Boundary.Fix{q}.SetName;
            
            attr = docNode.createAttribute('set'); %Create id attribute
            attr.setNodeValue(currentSetName); %Set id text
            parent_node.setAttributeNode(attr); %Add id attribute            
        end
        
        %Create a node set
        if isfield(FEB_struct.Boundary.Fix{q},'Set')
            currentSet=FEB_struct.Boundary.Fix{q}.Set;
            for q_node=1:1:numel(currentSet)
                node_node = docNode.createElement('node'); %create node entry
                node_node = parent_node.appendChild(node_node); %add node entry
                
                attr = docNode.createAttribute('id'); %Create id attribute
                attr.setNodeValue(sprintf('%u',currentSet(q_node))); %Set id text
                node_node.setAttributeNode(attr); %Add id attribute
            end
        end
    end
end

%% ADDING PRESCRIBED BC'S

if  isfield(FEB_struct.Boundary,'Prescribe')
    disp('----> Defining prescribe type boundary conditions');
    
    for q=1:1:numel(FEB_struct.Boundary.Prescribe) %For all the prescribe type BC's
        
        %Create prescribe section
        parent_node = docNode.createElement('prescribe');
        parent_node = boundaryNode.appendChild(parent_node);
        
        %Set bc attribute
        currentBC=FEB_struct.Boundary.Prescribe{q}.bc;
        attr = docNode.createAttribute('bc'); %Create id attribute
        attr.setNodeValue(currentBC); %Set id text
        parent_node.setAttributeNode(attr); %Add id attribute
        
        %Set optional scale parameter
        if isfield(FEB_struct.Boundary.Prescribe{q},'Scale')
            currentScale=FEB_struct.Boundary.Prescribe{q}.Scale;
            attr = docNode.createAttribute('scale'); %Create id attribute
            attr.setNodeValue(sprintf('%u',currentScale)); %Set id text
            parent_node.setAttributeNode(attr); %Add id attribute
        end        
        
        %Add set name attribute if present
        if isfield(FEB_struct.Boundary.Prescribe{q},'SetName')
            currentSetName=FEB_struct.Boundary.Prescribe{q}.SetName;
            
            attr = docNode.createAttribute('set'); %Create id attribute
            attr.setNodeValue(currentSetName); %Set id text
            parent_node.setAttributeNode(attr); %Add id attribute
            
        end
        
        %Add loadCurve attribute if present
        if isfield(FEB_struct.Boundary.Prescribe{q},'lc')
            currentLC=FEB_struct.Boundary.Prescribe{q}.lc;
            
            attr = docNode.createAttribute('lc'); %Create id attribute
            attr.setNodeValue(sprintf('%u',currentLC)); %Set id text
            parent_node.setAttributeNode(attr); %Add id attribute
            
        end
        
        %Add set name attribute if present
        if isfield(FEB_struct.Boundary.Prescribe{q},'Type')
            currentType=FEB_struct.Boundary.Prescribe{q}.Type;
            
            attr = docNode.createAttribute('type'); %Create id attribute
            attr.setNodeValue(currentType); %Set id text
            parent_node.setAttributeNode(attr); %Add id attribute
            
        end
        
        %Create node set
        if isfield(FEB_struct.Boundary.Prescribe{q},'Set')
            currentSet=FEB_struct.Boundary.Prescribe{q}.Set;
            for q_node=1:1:numel(currentSet)
                node_node = docNode.createElement('node'); %create node entry
                node_node = parent_node.appendChild(node_node); %add node entry
                
                attr = docNode.createAttribute('id'); %Create id attribute
                attr.setNodeValue(sprintf('%u',currentSet(q_node))); %Set id text
                node_node.setAttributeNode(attr); %Add id attribute
                
                currentValue=FEB_struct.Boundary.Prescribe{q}.nodeScale(q_node);
                node_node.appendChild(docNode.createTextNode(sprintf('%6.7e',currentValue))); %append data text child
            end
        end
    end
end


end