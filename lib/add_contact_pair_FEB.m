function docNode=add_contact_pair_FEB(docNode,FEB_struct,disp_opt)

% function docNode=add_contact_pair_FEB(docNode,FEB_struct,disp_opt)
% ------------------------------------------------------------------------
% This function adds FEBio contact definitions the input docNode defining
% an XML file suitable for generating .feb files.
%
%
% Change log;
% 2013/08/19 KMM: Changed .quad4 into .elements for contact surfaces specification.
%
% 2014/03/06 KMM: Added rigid_wall type contact and changed to swith statement
% for ContactType
%
% 2012/06/13, Kevin Mattheus Moerman
% kevinmoerman@hotmail.com

%%

%Get XML
febio_spec = docNode.getDocumentElement;

%Check for Boundary level
if docNode.getElementsByTagName('Boundary').getLength==0; %level does not exist yet
    ElementAddNode = docNode.createElement('Boundary');
    febio_spec.appendChild(ElementAddNode);
    BoundaryNode = docNode.getElementsByTagName('Boundary').item(0);
    IsBoundaryNodeNew=1;
else
    BoundaryNode = docNode.getElementsByTagName('Boundary').item(0);
    IsBoundaryNodeNew=0;
end


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
            if disp_opt==1;
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
                if disp_opt==1;
                    waitbar(i/n_steps);
                end
            end
            if disp_opt==1;
                close(hw);
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
                if disp_opt==1;
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
                    
                    i=i+1;
                    if disp_opt==1;
                        waitbar(i/n_steps);
                    end
                end
                if disp_opt==1;
                    close(hw);
                end
            end
    end
end
end
