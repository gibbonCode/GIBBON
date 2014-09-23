function docNode=add_constraint_section_FEB(docNode,FEB_struct)


% 29/08/2012, Kevin Mattheus Moerman
% kevinmoerman@gmail.com

%%

%Get XML
febio_spec = docNode.getDocumentElement;

%Check for Contraint level
if docNode.getElementsByTagName('Constraints').getLength==0; %level does not exist yet
    ElementAddNode = docNode.createElement('Constraints');
    febio_spec.appendChild(ElementAddNode);
    ControlNode = docNode.getElementsByTagName('Constraints').item(0);
    IsControlNodeNew=1;
else
    ControlNode = docNode.getElementsByTagName('Constraints').item(0);
    IsControlNodeNew=0;
end
docRootNode=ControlNode;
% 
% %Setting Control parameters
% controlProps=FEB_struct.Control.Properties;
% controlPropVals=FEB_struct.Control.Values; %default values
% for q=1:1:numel(controlProps)
%     prop_node = docNode.createElement(controlProps{q}); %create entry
%     prop_node = docRootNode.appendChild(prop_node); %add entry
%     t_form=repmat('%f, ',1,size(controlPropVals{q},2)); t_form=t_form(1:end-2);
%     prop_node.appendChild(docNode.createTextNode(sprintf(t_form,controlPropVals{q}))); %append data text child
% end
% 
% %Adding analysis type
% prop_node = docNode.createElement('analysis'); %create entry
% prop_node = docRootNode.appendChild(prop_node); %add entry
% attr = docNode.createAttribute('type'); %Create id attribute
% attr.setNodeValue(FEB_struct.Control.AnalysisType); %Set id text
% prop_node.setAttributeNode(attr); %Add id attribute
% 

%Rigid constraint entries
rigidId=FEB_struct.Constraints.RigidId;
prop_node = docNode.createElement('rigid_body'); %create entry
prop_node = docRootNode.appendChild(prop_node); %add entry
attr = docNode.createAttribute('mat'); %Create attribute
attr.setNodeValue(num2str(rigidId)); %Set text
prop_node.setAttributeNode(attr); %Add attribute

rigidProps=FEB_struct.Constraints.Properties;
rigidPropVals=FEB_struct.Constraints.Values;
rigidType=FEB_struct.Constraints.Type;
loadCurves=FEB_struct.Constraints.LoadCurveIds;

for q=1:1:numel(rigidProps)
    prop_prop_node = docNode.createElement(rigidProps{q}); %create entry
    prop_prop_node = prop_node.appendChild(prop_prop_node); %add entry
    
    attr = docNode.createAttribute('type'); %Create attribute
    attr.setNodeValue(rigidType); %Set text
    prop_prop_node.setAttributeNode(attr); %Add attribute
    
    attr = docNode.createAttribute('lc'); %Create attribute
    attr.setNodeValue(num2str(loadCurves(q))); %Set text
    prop_prop_node.setAttributeNode(attr); %Add attribute
    
    t_form=repmat('%f, ',1,size(rigidPropVals{q},2)); t_form=t_form(1:end-2);
    prop_prop_node.appendChild(docNode.createTextNode(sprintf(t_form,rigidPropVals{q}))); %append data text child
end

%% Add load curves
% 	<LoadData>
% 		<loadcurve id="1" type="smooth">
% 			<loadpoint>0,0</loadpoint>
% 			<loadpoint>1,1</loadpoint>
% 		</loadcurve>

%Check for LoadData level
if docNode.getElementsByTagName('LoadData').getLength==0; %level does not exist yet
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
% 
% for q=1:1:numel(rigidProps)
%     prop_prop_node = docNode.createElement(rigidProps{q}); %create entry
%     prop_prop_node = prop_node.appendChild(prop_prop_node); %add entry
%     
%     attr = docNode.createAttribute('type'); %Create attribute
%     attr.setNodeValue(rigidType); %Set text
%     prop_prop_node.setAttributeNode(attr); %Add attribute
%     
%     attr = docNode.createAttribute('lc'); %Create attribute
%     attr.setNodeValue(num2str(loadCurves(q))); %Set text
%     prop_prop_node.setAttributeNode(attr); %Add attribute
%     
%     t_form=repmat('%f, ',1,size(rigidPropVals{q},2)); t_form=t_form(1:end-2);
%     prop_prop_node.appendChild(docNode.createTextNode(sprintf(t_form,rigidPropVals{q}))); %append data text child
% end

end
% write_XML_no_extra_lines(savename,docNode)% Saving XML file