function [docNode]=addLoadsLevel_FEB(docNode,FEB_struct)


%%

disp('Adding Loads level');

rootNode = docNode.getDocumentElement;

loadsNode = docNode.createElement('Loads');
loadsNode = rootNode.appendChild(loadsNode);

%% Adding surface loads
if isfield(FEB_struct.Loads,'Surface_load')
    for q_sl=1:1:numel(FEB_struct.Loads.Surface_load) %For all surface_loads
        
        disp('----> Defining surface loads')
        
        %Create surface_load level
        surfaceLoadNode = docNode.createElement('surface_load');
        surfaceLoadNode = loadsNode.appendChild(surfaceLoadNode);
        
        %Set load type attribute
        currentLoadType=FEB_struct.Loads.Surface_load{q_sl}.Type;
        attr = docNode.createAttribute('type'); %Create id attribute
        attr.setNodeValue(currentLoadType); %Set attribute text
        surfaceLoadNode.setAttributeNode(attr); %Add id attribute
        
        %Create lcPar level
        currentLcParType=FEB_struct.Loads.Surface_load{q_sl}.lcPar;
        lcParNode = docNode.createElement(currentLcParType);
        lcParNode = surfaceLoadNode.appendChild(lcParNode);
        
        %Set lcPar value (sometimes scale)
        lcParNode.appendChild(docNode.createTextNode(sprintf('%10.6e',FEB_struct.Loads.Surface_load{q_sl}.lcParValue))); %append data text child
        
        %Set lcPar lc attribute
        if isfield(FEB_struct.Loads.Surface_load{q_sl},'lc')
            attr = docNode.createAttribute('lc'); %Create id attribute
            attr.setNodeValue(sprintf('%u',FEB_struct.Loads.Surface_load{q_sl}.lc)); %Set attribute text
            lcParNode.setAttributeNode(attr); %Add id attribute
        end
        
        %Set other properties if present
        if isfield(FEB_struct.Loads.Surface_load{q_sl},'Properties')
            surfaceLoadProperties=FEB_struct.Loads.Surface_load{q_sl}.Properties;
            surfaceLoadValues=FEB_struct.Loads.Surface_load{q_sl}.Values;
            
            for q=1:1:numel(surfaceLoadProperties)
                prop_node = docNode.createElement(surfaceLoadProperties{q}); %create entry
                prop_node = surfaceLoadNode.appendChild(prop_node); %add entry
                t_form=repmat('%10.6e, ',1,size(surfaceLoadValues{q},2)); t_form=t_form(1:end-2);
                prop_node.appendChild(docNode.createTextNode(sprintf(t_form,surfaceLoadValues{q}))); %append data text child
            end
        end
        
        %Create surface_load level
        disp('----> Adding loaded surface');
        surfaceNode = docNode.createElement('surface');
        surfaceNode = surfaceLoadNode.appendChild(surfaceNode);
        if isfield(FEB_struct.Loads.Surface_load{q_sl},'SetName') && ~isfield(FEB_struct.Loads.Surface_load{q_sl},'Set')
            %Define set name attribute
            attr = docNode.createAttribute('set'); %Create id attribute
            attr.setNodeValue(FEB_struct.Loads.Surface_load{q_sl}.SetName); %Set attribute text
            surfaceNode.setAttributeNode(attr); %Add id attribute
        elseif isfield(FEB_struct.Loads.Surface_load{q_sl},'Set') && ~isfield(FEB_struct.Loads.Surface_load{q_sl},'SetName')
            
            %Specifying the surface elements
            E=FEB_struct.Loads.Surface_load{q_sl}.Set; %The set of surface elements
            n_steps=size(E,1); %Number of surface elements
            
            %Define waitbar
            if FEB_struct.disp_opt==1;
                hw = waitbar(0,'Adding contact element entries....');
            end
            for q_e=1:1:n_steps
                
                %Set surface type based on number of nodes
                switch size(E,2)
                    case 4 %quad4
                        E_type='quad4';
                    case 3 %tri3
                        E_type='tri3';
                    case 6 %tri6
                        E_type='tri6';
                end
                
                %Create surface type entry
                element_node = docNode.createElement(E_type); %create element entry
                element_node = surfaceNode.appendChild(element_node); %add element entry
                
                %Set id
                attr = docNode.createAttribute('id'); %Create id attribute
                attr.setNodeValue(num2str(q_e)); %Set id text
                element_node.setAttributeNode(attr); %Add id attribute
                t_form=repmat('%u, ',1,size(E,2)); t_form=t_form(1:end-2); %Text form
                element_node.appendChild(docNode.createTextNode(sprintf(t_form,E(q_e,:)))); %append data text child
                
                %Increment waitbar
                if FEB_struct.disp_opt==1 && rem(q_e,round(n_steps/10))==0;
                    waitbar(q_e/n_steps);
                end
            end
            
            if FEB_struct.disp_opt==1;
                close(hw); drawnow;
            end
            
        else
            error('Specify either SetName or Set not both');
        end
    end
end

%% Adding nodal loads

if  isfield(FEB_struct.Loads,'Nodal_load')
    
    disp('----> Defining node loads')
    
    for q1=1:1:numel(FEB_struct.Loads.Nodal_load) %For all Nodal_loads
        
        %Create nodal_load level
        nodal_load_node = docNode.createElement('nodal_load'); %create entry
        nodal_load_node = loadsNode.appendChild(nodal_load_node); %add entry
        
        %Add bc attribute
        currentBC=FEB_struct.Loads.Nodal_load{q1}.bc;
        attr = docNode.createAttribute('bc'); %Create id attribute
        attr.setNodeValue(currentBC); %Set id text
        nodal_load_node.setAttributeNode(attr); %Add id attribute
        
        %Set lc attribute
        if isfield(FEB_struct.Loads.Nodal_load{q1},'lc')
            attr = docNode.createAttribute('lc'); %Create id attribute
            attr.setNodeValue(sprintf('%u',FEB_struct.Loads.Nodal_load{q1}.lc)); %Set attribute text
            nodal_load_node.setAttributeNode(attr); %Add id attribute
        end
        
        %Create node entries
        nodeIdSet=FEB_struct.Loads.Nodal_load{q1}.Set;
        nodeValues=FEB_struct.Loads.Nodal_load{q1}.nodeScale;
        
        if FEB_struct.disp_opt==1;
            hw = waitbar(0,'Adding nodal loads....');
        end
        n_steps=size(nodeIdSet,1);
        
        disp('----> Adding nodal loads');
        
        for q2=1:1:n_steps
            
            node_node = docNode.createElement('node'); %create node entry
            node_node = nodal_load_node.appendChild(node_node); %add node entry
            
            attr = docNode.createAttribute('id'); %Create id attribute
            attr.setNodeValue(sprintf('%u',nodeIdSet(q2))); %Set id text
            node_node.setAttributeNode(attr); %Add id attribute
            
            node_node.appendChild(docNode.createTextNode(sprintf('%6.7e',nodeValues(q2)))); %append data text child
            
            if FEB_struct.disp_opt==1 && rem(q2,round(n_steps/10))==0;
                waitbar(q2/n_steps);
            end
        end
        if FEB_struct.disp_opt==1;
            close(hw); drawnow;
        end
        
    end
end


