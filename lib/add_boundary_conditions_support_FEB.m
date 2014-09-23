function docNode=add_boundary_conditions_support_FEB(docNode,FEB_struct)

% %function docNode=add_boundary_conditions_support_FEB(docNode,FEB_struct)
%
%
%
% CHANGE LOG:
% 29/08/2012, Kevin Mattheus Moerman
% 02/12/2013, Kevin Mattheus Moerman, added prescribed displacement boundary conditions
% 05/15/2014, Kevin Mattheus Moerman, added pressure load boundary conditions
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% ------------------------------------------------------------------------

%%

disp('Adding Boundary level')

%Get XML
febio_spec = docNode.getDocumentElement;

%Check for Boundary level
% if docNode.getElementsByTagName('Boundary').getLength==0; %level does not exist yet
%     ElementAddNode = docNode.createElement('Boundary');
%     febio_spec.appendChild(ElementAddNode);
%     BoundaryNode = docNode.getElementsByTagName('Boundary').item(0);
%     IsBoundaryNodeNew=1;
% else
%     BoundaryNode = docNode.getElementsByTagName('Boundary').item(0);
%     IsBoundaryNodeNew=0;
% end

if  isfield(FEB_struct,'Boundary')
    
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
    
    %% ADDING FIX BC's
    docRootNode=BoundaryNode;
    
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
        if BoundaryNode.getElementsByTagName('prescribe').getLength==0; %level does not exist yet
            parent_node = docNode.createElement('prescribe');
            parent_node = docRootNode.appendChild(parent_node);
            IsContactNodeNew=1;
        else
            IsContactNodeNew=0; %TO DO not treating this case yet
        end
        docRootNode =  BoundaryNode.getElementsByTagName('prescribe').item(0);
        
        loadCurves=FEB_struct.Boundary.LoadCurveIds;
        prePropVals=FEB_struct.Boundary.PrescribeValues;
        
        for q1=1:1:numel(FEB_struct.Boundary.PrescribeList);
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
                
                t_form=repmat('%f, ',1,size(prePropValSet(q2),2));
                t_form=t_form(1:end-2);
                node_node.appendChild(docNode.createTextNode(sprintf(t_form,prePropValSet(q2)))); %append data text child
            end
        end
        
    end
end

%%
disp('----> Defining loads');

if  isfield(FEB_struct,'Loads')
    
    %     			<tri3 id="3097" lc="1" scale="0.018">  4480,  4067,  4718</tri3>
    % 		</pressure>
    
    %Check for Loads level
    if docNode.getElementsByTagName('<Loads>').getLength==0; %level does not exist yet
        ElementAddNode = docNode.createElement('Loads');
        febio_spec.appendChild(ElementAddNode);
        ControlNode = docNode.getElementsByTagName('Loads').item(0);
        IsControlNodeNew=1;
    else
        ControlNode = docNode.getElementsByTagName('Loads').item(0);
        IsControlNodeNew=0;
    end
    docRootNode=ControlNode;
    
    if  isfield(FEB_struct.Loads,'Pressure')
        
        for qPres=1:1:numel(FEB_struct.Loads.Pressure.Elements)
            
            prop_node = docNode.createElement('pressure'); %create entry
            prop_node = docRootNode.appendChild(prop_node); %add entry
            
            disp('----> Adding pressure surface');
            E_types={'quad4','tri3','tri6'};
            
            E=FEB_struct.Loads.Pressure.Elements{qPres};
            
            hw = waitbar(0,'Adding pressure surface element entries....');
            
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
                element_node = prop_node.appendChild(element_node); %add element entry
                
                attr = docNode.createAttribute('id'); %Create id attribute
                attr.setNodeValue(num2str(i)); %Set id text
                element_node.setAttributeNode(attr); %Add id attribute
                
                attr = docNode.createAttribute('lc'); %Create lc attribute
                attr.setNodeValue(num2str(FEB_struct.Loads.Pressure.LoadCurveIds(qPres))); %Set id text
                element_node.setAttributeNode(attr); %Add id attribute
                
                attr = docNode.createAttribute('scale'); %Create lc attribute
                attr.setNodeValue(sprintf('%f',FEB_struct.Loads.Pressure.Scale(qPres))); %Set id text
                element_node.setAttributeNode(attr); %Add id attribute
                
                t_form=repmat('%d, ',1,size(E,2)); t_form=t_form(1:end-2);
                element_node.appendChild(docNode.createTextNode(sprintf(t_form,E(i,:)))); %append data text child
                
                waitbar(i/n_steps);
                
            end            
            close(hw);           
        end
    end
end

%% Add load curves
disp('----> Defining load curves');

if  isfield(FEB_struct,'LoadData')
    % 	<LoadData>
    % 		<loadcurve id="1" type="smooth">
    % 			<loadpoint>0,0</loadpoint>
    % 			<loadpoint>1,1</loadpoint>
    % 		</loadcurve>
    
    %Check for LoadData level
    if docNode.getElementsByTagName('<LoadData>').getLength==0; %level does not exist yet
        ElementAddNode = docNode.createElement('LoadData');
        febio_spec.appendChild(ElementAddNode);
        ControlNode = docNode.getElementsByTagName('LoadData').item(0);
        IsControlNodeNew=1;
    else
        ControlNode = docNode.getElementsByTagName('LoadData').item(0);
        IsControlNodeNew=0;
    end
    docRootNode=ControlNode;
    
    for q=1:1:numel(FEB_struct.LoadData.LoadCurves.id)
        
        loadPoints=FEB_struct.LoadData.LoadCurves.loadPoints{q};
        
        prop_node = docNode.createElement('loadcurve'); %create entry
        prop_node = docRootNode.appendChild(prop_node); %add entry
        
        attr = docNode.createAttribute('id'); %Create attribute
        attr.setNodeValue(num2str(FEB_struct.LoadData.LoadCurves.id(q))); %Set text
        prop_node.setAttributeNode(attr); %Add attribute
        
        attr = docNode.createAttribute('type'); %Create attribute
        attr.setNodeValue(FEB_struct.LoadData.LoadCurves.type{q}); %Set text
        prop_node.setAttributeNode(attr); %Add attribute
        
        for q2=1:1:size(loadPoints,1)
            loadPoint=loadPoints(q2,:);
            prop_prop_node = docNode.createElement('loadpoint'); %create entry
            prop_prop_node = prop_node.appendChild(prop_prop_node); %add entry
            
            t_form=repmat('%f, ',1,size(loadPoint,2)); t_form=t_form(1:end-2);
            prop_prop_node.appendChild(docNode.createTextNode(sprintf(t_form,loadPoint))); %append data text child
        end
    end
end

end

