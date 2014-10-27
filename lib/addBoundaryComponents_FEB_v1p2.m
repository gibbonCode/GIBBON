function [docNode]=addBoundaryComponents_FEB_v1p2(docNode,BoundaryNode,FEB_struct)

docRootNode=BoundaryNode;

%% ADDING FIX BC's

if  isfield(FEB_struct.Boundary,'FixList')
    disp('----> Defining fix type boundary conditions');
    
    %Check for fix level
    if BoundaryNode.getElementsByTagName('fix').getLength==0; %level does not exist yet
        parent_node = docNode.createElement('fix');
        parent_node = docRootNode.appendChild(parent_node);
        IsContactNodeNew=1;
    else
        IsContactNodeNew=0; %TO DO not treating this case yet
    end
    docRootNode =  BoundaryNode.getElementsByTagName('fix').item(0);
    
    for q1=1:1:numel(FEB_struct.Boundary.FixList);
        nodeList=FEB_struct.Boundary.FixList{q1};
        fixType=FEB_struct.Boundary.FixType{q1};
        for q2=1:1:numel(nodeList);
            node_node = docNode.createElement('node'); %create node entry
            node_node = parent_node.appendChild(node_node); %add node entry
            attr = docNode.createAttribute('id'); %Create id attribute
            attr.setNodeValue(num2str(nodeList(q2))); %Set id text
            node_node.setAttributeNode(attr); %Add id attribute
            attr = docNode.createAttribute('bc'); %Create id attribute
            attr.setNodeValue(fixType); %Set id text
            node_node.setAttributeNode(attr); %Add id attribute
        end
    end
    
end

%% ADDING PRESCRIBED BC'S

if  isfield(FEB_struct.Boundary,'PrescribeList')
    disp('----> Defining prescribe type boundary conditions');
    
    %%
    docRootNode=BoundaryNode;
    %Check for prescribe level
    %         if BoundaryNode.getElementsByTagName('prescribe').getLength==0; %level does not exist yet
    %             parent_node = docNode.createElement('prescribe');
    %             parent_node = docRootNode.appendChild(parent_node);
    %             IsContactNodeNew=1;
    %         else
    %             IsContactNodeNew=0; %TO DO not treating this case yet
    %         end
    %         docRootNode =  BoundaryNode.getElementsByTagName('prescribe').item(0);
    
    loadCurves=FEB_struct.Boundary.LoadCurveIds;
    prePropVals=FEB_struct.Boundary.PrescribeValues;
    
    for q1=1:1:numel(FEB_struct.Boundary.PrescribeList);
        
        %Create prescribe section
        parent_node = docNode.createElement('prescribe');
        parent_node = docRootNode.appendChild(parent_node);
        
        %Define prescribe type if specified
        if isfield(FEB_struct.Boundary,'PrescribeTypes')
            attr = docNode.createAttribute('type'); %Create id attribute
            attr.setNodeValue(FEB_struct.Boundary.PrescribeTypes(q1)); %Set id text
            parent_node.setAttributeNode(attr); %Add id attribute
        end
        nodeList=FEB_struct.Boundary.PrescribeList{q1};
        preType=FEB_struct.Boundary.PrescribeType{q1};
        prePropValSet=prePropVals{q1};
        for q2=1:1:numel(nodeList);
            
            node_node = docNode.createElement('node'); %create node entry
            node_node = parent_node.appendChild(node_node); %add node entry
            
            attr = docNode.createAttribute('id'); %Create id attribute
            attr.setNodeValue(num2str(nodeList(q2))); %Set id text
            node_node.setAttributeNode(attr); %Add id attribute
            
            attr = docNode.createAttribute('bc'); %Create id attribute
            attr.setNodeValue(preType); %Set id text
            node_node.setAttributeNode(attr); %Add id attribute
            
            attr = docNode.createAttribute('lc'); %Create attribute
            attr.setNodeValue(num2str(loadCurves(q1))); %Set text
            node_node.setAttributeNode(attr); %Add attribute
            
            t_form=repmat('%f, ',1,size(prePropValSet(q2,:),2));
            t_form=t_form(1:end-2);                       
      
            node_node.appendChild(docNode.createTextNode(sprintf(t_form,prePropValSet(q2,:)))); %append data text child
        end
    end
    
end

%% ADDING CONTACT

if  isfield(FEB_struct.Boundary,'Contact')
    
    disp('----> Defining contact')
    
    for q_set=1:1:numel(FEB_struct.Boundary.Contact)
        docRootNode=BoundaryNode;
        
        %Add (another if already present) contact level
        parent_node = docNode.createElement('contact');
        parent_node = docRootNode.appendChild(parent_node);
        numContactSets=BoundaryNode.getElementsByTagName('contact').getLength; %Current number of contact sets
        docRootNode =  BoundaryNode.getElementsByTagName('contact').item(numContactSets-1);
        
        %Setting Contact Type
        ContactType=FEB_struct.Boundary.Contact{q_set}.Type;
        attr = docNode.createAttribute('type'); %Create attribute
        attr.setNodeValue(ContactType); %Set text
        docRootNode.setAttributeNode(attr); %Add attribute
        
        switch ContactType
            case 'rigid'
                disp('----> Adding rigid body nodes for rigid contact');
                rigidNodeList=FEB_struct.Boundary.Contact{q_set}.RigidNodeList;
                rigidIdList=FEB_struct.Boundary.Contact{q_set}.RigidIdList;
                for q=1:1:numel(rigidNodeList)
                    element_node = docNode.createElement('node'); %create element entry
                    element_node = parent_node.appendChild(element_node); %add element entry
                    
                    attr = docNode.createAttribute('id'); %Create id attribute
                    attr.setNodeValue(num2str(rigidNodeList(q))); %Set id text
                    element_node.setAttributeNode(attr); %Add id attribute
                    
                    attr = docNode.createAttribute('rb'); %Create id attribute
                    attr.setNodeValue(num2str(rigidIdList(q))); %Set id text
                    element_node.setAttributeNode(attr); %Add id attribute
                end
            case 'rigid_wall'
                %Setting Contact parameters
                disp('----> Setting contact parameters');
                contactProperties=FEB_struct.Boundary.Contact{q_set}.Properties;
                contactValues=FEB_struct.Boundary.Contact{q_set}.Values;
                for q=1:1:numel(contactProperties)
                    prop_node = docNode.createElement(contactProperties{q}); %create entry
                    prop_node = docRootNode.appendChild(prop_node); %add entry
                    t_form=repmat('%10.6e, ',1,size(contactValues{q},2)); t_form=t_form(1:end-2);
                    prop_node.appendChild(docNode.createTextNode(sprintf(t_form,contactValues{q}))); %append data text child
                end
                
                %Defining rigid wall plane
                disp('----> Defining rigid wall plane');
                prop_node = docNode.createElement('plane'); %create entry
                prop_node = docRootNode.appendChild(prop_node); %add entry
                t_form=repmat('%10.6e, ',1,size(FEB_struct.Boundary.Contact{q_set}.planePar{1},2)); t_form=t_form(1:end-2);
                prop_node.appendChild(docNode.createTextNode(sprintf(t_form,FEB_struct.Boundary.Contact{q_set}.planePar{1}))); %append data text child
                
                attr = docNode.createAttribute('lc'); %Create attribute
                attr.setNodeValue(num2str(FEB_struct.Boundary.Contact{q_set}.loadCurve)); %Set text
                prop_node.setAttributeNode(attr); %Add attribute
                
                %% Adding contact surface
                
                disp('----> Adding contact surface');
                E_types={'quad4','tri3','tri6'};
                %         docRootNode =  BoundaryNode.getElementsByTagName('contact').item(numContactSets-1);
                
                parent_node = docNode.createElement('surface');
                parent_node = docRootNode.appendChild(parent_node);
                %         docRootNode = BoundaryNode.getElementsByTagName('surface').item(0);
                
                E=FEB_struct.Boundary.Contact{q_set}.Surfaces.elements{1};
                if FEB_struct.disp_opt==1;
                    hw = waitbar(0,'Adding contact element entries....');
                end
                n_steps=size(E,1);
                
                for i=1:1:n_steps
                    switch size(E,2)
                        case 4 %quad4
                            E_type=E_types{1};
                        case 3 %tri3
                            E_type=E_types{2};
                        case 6 %tri6
                            E_type=E_types{3};
                    end
                    
                    element_node = docNode.createElement(E_type); %create element entry
                    element_node = parent_node.appendChild(element_node); %add element entry
                    
                    attr = docNode.createAttribute('id'); %Create id attribute
                    attr.setNodeValue(num2str(i)); %Set id text
                    element_node.setAttributeNode(attr); %Add id attribute
                    
                    t_form=repmat('%d, ',1,size(E,2)); t_form=t_form(1:end-2);
                    element_node.appendChild(docNode.createTextNode(sprintf(t_form,E(i,:)))); %append data text child
                    
                    %                 i=i+1;
                    if FEB_struct.disp_opt==1 && rem(i,round(n_steps/10))==0;
                        waitbar(i/n_steps);
                    end
                end
                if FEB_struct.disp_opt==1;
                    close(hw); drawnow;
                end
                
            otherwise %Sliding, Tied or Sticky assumed
                disp('----> Setting contact parameters');
                %Setting Contact parameters
                contactProperties=FEB_struct.Boundary.Contact{q_set}.Properties;
                contactValues=FEB_struct.Boundary.Contact{q_set}.Values;
                
                for q=1:1:numel(contactProperties)
                    prop_node = docNode.createElement(contactProperties{q}); %create entry
                    prop_node = docRootNode.appendChild(prop_node); %add entry
                    t_form=repmat('%10.6e, ',1,size(contactValues{q},2)); t_form=t_form(1:end-2);
                    prop_node.appendChild(docNode.createTextNode(sprintf(t_form,contactValues{q}))); %append data text child
                end
                
                %% Adding contact surface pairs
                E_types={'quad4','tri3','tri6'};
                for q=1:1:2
                    docRootNode =  BoundaryNode.getElementsByTagName('contact').item(numContactSets-1);
                    SurfaceType=FEB_struct.Boundary.Contact{q_set}.Surfaces.Type{q};
                    parent_node = docNode.createElement('surface');
                    parent_node = docRootNode.appendChild(parent_node);
                    attr = docNode.createAttribute('type'); %Create attribute
                    attr.setNodeValue(SurfaceType); %Set text
                    parent_node.setAttributeNode(attr); %Add attribute
                    
                    E=FEB_struct.Boundary.Contact{q_set}.Surfaces.elements{q};
                    if FEB_struct.disp_opt==1;
                        hw = waitbar(0,'Adding contact element entries....');
                    end
                    n_steps=size(E,1);
                    if q==1
                        disp('----> Adding MASTER contact element entries')
                    else
                        disp('----> Adding SLAVE contact element entries')
                    end
                    
                    for i=1:1:n_steps
                        switch size(E,2)
                            case 4 %quad4
                                E_type=E_types{1};
                            case 3 %tri3
                                E_type=E_types{2};
                            case 6 %tri6
                                E_type=E_types{3};
                        end
                        
                        element_node = docNode.createElement(E_type); %create element entry
                        element_node = parent_node.appendChild(element_node); %add element entry
                        
                        attr = docNode.createAttribute('id'); %Create id attribute
                        attr.setNodeValue(num2str(i)); %Set id text
                        element_node.setAttributeNode(attr); %Add id attribute
                        
                        t_form=repmat('%d, ',1,size(E,2)); t_form=t_form(1:end-2);
                        element_node.appendChild(docNode.createTextNode(sprintf(t_form,E(i,:)))); %append data text child
                        
                        %                             i=i+1;
                        if FEB_struct.disp_opt==1 && rem(i,round(n_steps/10))==0;
                            waitbar(i/n_steps);
                        end
                    end
                    if FEB_struct.disp_opt==1;
                        close(hw); drawnow;
                    end
                end
        end
    end
end