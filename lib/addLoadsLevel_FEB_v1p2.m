function [docNode]=addLoadsLevel_FEB_v1p2(docNode,FEB_struct)


%%

disp('Adding Loads level')

febio_spec = docNode.getDocumentElement;

if  isfield(FEB_struct,'Loads')
    
    %Check for Loads level
    if docNode.getElementsByTagName('<Loads>').getLength==0; %level does not exist yet
        ElementAddNode = docNode.createElement('Loads');
        febio_spec.appendChild(ElementAddNode);
        ControlNode = docNode.getElementsByTagName('Loads').item(0);
%         IsControlNodeNew=1;
    else
        ControlNode = docNode.getElementsByTagName('Loads').item(0);
%         IsControlNodeNew=0;
    end
    docRootNode=ControlNode;
    
    %% Adding pressure load
    if  isfield(FEB_struct.Loads,'Pressure')
        disp('----> Pressure load')
        for q1=1:1:numel(FEB_struct.Loads.Pressure.Elements)
            
            prop_node = docNode.createElement('pressure'); %create entry
            prop_node = docRootNode.appendChild(prop_node); %add entry
            
            disp('----> Adding pressure surface');
            E_types={'quad4','tri3','tri6'};
            
            E=FEB_struct.Loads.Pressure.Elements{q1}; %Element set
            P=FEB_struct.Loads.Pressure.Scale{q1}; %Pressure set
            
            if FEB_struct.disp_opt==1;
                hw = waitbar(0,'Adding pressure surface element entries....');
            end
            n_steps=size(E,1);
            
            for q2=1:1:n_steps
                switch size(E,2)
                    case 4 %quad4
                        E_type=E_types{1};
                    case 3 %tri3
                        E_type=E_types{2};
                    case 6 %tri6
                        E_type=E_types{3};
                end
                
                element_node = docNode.createElement(E_type); %create element entry
                element_node = prop_node.appendChild(element_node); %add element entry
                
                attr = docNode.createAttribute('id'); %Create id attribute
                attr.setNodeValue(num2str(q2)); %Set id text
                element_node.setAttributeNode(attr); %Add id attribute
                
                attr = docNode.createAttribute('lc'); %Create lc attribute
                attr.setNodeValue(num2str(FEB_struct.Loads.Pressure.LoadCurveIds(q1))); %Set id text
                element_node.setAttributeNode(attr); %Add id attribute
                
                attr = docNode.createAttribute('scale'); %Create lc attribute
                attr.setNodeValue(sprintf('%f',P(q2))); %Set id text
                element_node.setAttributeNode(attr); %Add id attribute
                
                t_form=repmat('%d, ',1,size(E,2)); t_form=t_form(1:end-2);
                element_node.appendChild(docNode.createTextNode(sprintf(t_form,E(q2,:)))); %append data text child
                
                if FEB_struct.disp_opt==1 && rem(q2,round(n_steps/10))==0;                   
                    waitbar(q2/n_steps);
                end
            end
            if FEB_struct.disp_opt==1;                
                close(hw); drawnow; 
            end
        end
    end
    
    %% Adding force loads
    if  isfield(FEB_struct.Loads,'Force')
        disp('----> Force load')
        for q1=1:1:numel(FEB_struct.Loads.Force.Nodes)
            
            prop_node = docNode.createElement('force'); %create entry
            prop_node = docRootNode.appendChild(prop_node); %add entry
            
            disp('----> Adding nodal forces');
            
            %%
            
            N=FEB_struct.Loads.Force.Nodes{q1};
            F=FEB_struct.Loads.Force.PrescribeValues{q1};
            
            if FEB_struct.disp_opt==1;
                hw = waitbar(0,'Adding nodal forces....');
            end
            n_steps=size(N,1);
            
            for q2=1:1:n_steps
                
                node_node = docNode.createElement('node'); %create node entry
                node_node = prop_node.appendChild(node_node); %add node entry
                
                attr = docNode.createAttribute('id'); %Create id attribute
                attr.setNodeValue(num2str(N(q2))); %Set id text
                node_node.setAttributeNode(attr); %Add id attribute
                
                attr = docNode.createAttribute('bc'); %Create bc attribute
                attr.setNodeValue(FEB_struct.Loads.Force.PrescribeType(q1)); %Set id text
                node_node.setAttributeNode(attr); %Add id attribute
                
                attr = docNode.createAttribute('lc'); %Create lc attribute
                attr.setNodeValue(num2str(FEB_struct.Loads.Force.LoadCurveIds(q1))); %Set id text
                node_node.setAttributeNode(attr); %Add id attribute
                
                node_node.appendChild(docNode.createTextNode(sprintf('%6.7e',F(q2)))); %append data text child
                
                if FEB_struct.disp_opt==1 && rem(q2,round(n_steps/10))==0;                   
                    waitbar(q2/n_steps);
                end
            end
            if FEB_struct.disp_opt==1;                
                close(hw); drawnow; 
            end

        end
        
    end
    
end
