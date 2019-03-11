function [docNode]=addContactComponents_FEB(docNode,contactNode,FEB_struct)

%%

if ~isfield(FEB_struct,'disp_opt');
    FEB_struct.disp_opt=0; 
end

disp('----> Defining contact')

for q_contact=1:1:numel(FEB_struct.Contact)
    
    %Create contact field within Contact
    contactContactNode = docNode.createElement('contact');
    contactContactNode = contactNode.appendChild(contactContactNode);
    
    %Setting Contact Type
    contactType=FEB_struct.Contact{q_contact}.Type;
    attr = docNode.createAttribute('type'); %Create attribute
    attr.setNodeValue(contactType); %Set text
    contactContactNode.setAttributeNode(attr); %Add attribute
    
    switch contactType
        case 'rigid'
            disp('----> Adding rigid body nodes for rigid contact');
            
            rigidNodeList=FEB_struct.Contact{q_contact}.RigidNodeList;
            rigidIdList=FEB_struct.Contact{q_contact}.RigidIdList;
            
            for q=1:1:numel(rigidNodeList)
                element_node = docNode.createElement('node'); %create element entry
                element_node = contactContactNode.appendChild(element_node); %add element entry
                
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
            parameterStruct=FEB_struct.Contact{q_contact};
            [docNode]=setContactProperties(docNode,contactContactNode,parameterStruct);
            
            %Defining rigid wall plane
            disp('----> Defining rigid wall plane');
            prop_node = docNode.createElement('plane'); %create entry
            prop_node = contactContactNode.appendChild(prop_node); %add entry
            t_form=repmat('%10.6e, ',1,size(FEB_struct.Contact{q_contact}.planePar{1},2)); t_form=t_form(1:end-2);
            prop_node.appendChild(docNode.createTextNode(sprintf(t_form,FEB_struct.Contact{q_contact}.planePar{1}))); %append data text child
            
            attr = docNode.createAttribute('lc'); %Create attribute
            attr.setNodeValue(num2str(FEB_struct.Contact{q_contact}.loadCurve)); %Set text
            prop_node.setAttributeNode(attr); %Add attribute
            
            % Adding contact surface
            disp('----> Defining contact surface');                        
            subSurfStruct=FEB_struct.Contact{q_contact}.Surface;
            subSurfStruct.disp_opt=FEB_struct.disp_opt;
            [docNode]=addSurfaceSpec(docNode,contactContactNode,subSurfStruct);
            
        otherwise %Sliding, Tied or Sticky assumed                        
            %Setting Contact parameters
            disp('----> Setting contact parameters');
            parameterStruct=FEB_struct.Contact{q_contact};
            [docNode]=setContactProperties(docNode,contactContactNode,parameterStruct);
            
            % Adding contact surface pairs
            disp('----> Defining contact surface pair');
            
            if numel(FEB_struct.Contact{q_contact}.Surface)~=2
                error('Wrong number of surface sets provided. Two sets are required.');
            end
            
            for q=1:1:2
                subSurfStruct=FEB_struct.Contact{q_contact}.Surface{q};
                subSurfStruct.disp_opt=FEB_struct.disp_opt;
                [docNode]=addSurfaceSpec(docNode,contactContactNode,subSurfStruct);
            end
    end
end

end

function [docNode]=setContactProperties(docNode,contactContactNode,parameterStruct)

contactProperties=parameterStruct.Properties;
contactValues=parameterStruct.Values;

for q=1:1:numel(contactProperties)
    prop_node = docNode.createElement(contactProperties{q}); %create entry
    prop_node = contactContactNode.appendChild(prop_node); %add entry
    t_form=repmat('%10.6e, ',1,size(contactValues{q},2)); t_form=t_form(1:end-2);
    prop_node.appendChild(docNode.createTextNode(sprintf(t_form,contactValues{q}))); %append data text child
end

end

function [docNode]=addSurfaceSpec(docNode,contactContactNode,subSurfStruct)

%Add surface level
surface_node = docNode.createElement('surface');
surface_node = contactContactNode.appendChild(surface_node);

%Add surface type (master or slave)
currentSurfaceType=subSurfStruct.Type;

attr = docNode.createAttribute('type'); %Create attribute
attr.setNodeValue(currentSurfaceType); %Set text
surface_node.setAttributeNode(attr); %Add attribute

%Specify surfaces using set or setName
if isfield(subSurfStruct,'Set') && ~isfield(subSurfStruct,'SetName') %A set is specified
    
    %Specifying the surface elements
    E=subSurfStruct.Set; %The set of surface elements    
    n_steps=size(E,1); %Number of surface elements
    
    %Define waitbar
    if subSurfStruct.disp_opt==1;
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
        element_node = surface_node.appendChild(element_node); %add element entry
        
        %Set id
        attr = docNode.createAttribute('id'); %Create id attribute
        attr.setNodeValue(num2str(q_e)); %Set id text
        element_node.setAttributeNode(attr); %Add id attribute        
        t_form=repmat('%u, ',1,size(E,2)); t_form=t_form(1:end-2); %Text form
        element_node.appendChild(docNode.createTextNode(sprintf(t_form,E(q_e,:)))); %append data text child
        
        %Increment waitbar
        if subSurfStruct.disp_opt==1 && rem(q_e,round(n_steps/10))==0;
            waitbar(q_e/n_steps);
        end
    end
    
    if subSurfStruct.disp_opt==1;
        close(hw); drawnow;
    end
    
elseif isfield(subSurfStruct,'SetName') && ~isfield(subSurfStruct,'Set') %A set name is given referring to the Surface section
    %Specifying the name of the set
    currentSurfaceSetName=subSurfStruct.SetName;
    attr = docNode.createAttribute('set'); %Create attribute
    attr.setNodeValue(currentSurfaceSetName); %Set text
    surface_node.setAttributeNode(attr); %Add attribute
else
    error('Cannot specify both a Set and a SetName for a contact surface');
end


end

 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
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
