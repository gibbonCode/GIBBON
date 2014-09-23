function docNode=add_control_section_FEB(docNode,FEB_struct)


%         
%
% 14/06/2012, Kevin Mattheus Moerman
% kevinmoerman@gmail.com

%%

%Get XML
febio_spec = docNode.getDocumentElement;

%Check for Control level
if docNode.getElementsByTagName('Control').getLength==0; %level does not exist yet
    ElementAddNode = docNode.createElement('Control');
    febio_spec.appendChild(ElementAddNode);
    ControlNode = docNode.getElementsByTagName('Control').item(0);
    IsControlNodeNew=1;
else
    ControlNode = docNode.getElementsByTagName('Control').item(0);
    IsControlNodeNew=0;
end
docRootNode=ControlNode;

%Setting Control parameters
controlProps=FEB_struct.Control.Properties;
controlPropVals=FEB_struct.Control.Values; %default values
for q=1:1:numel(controlProps)
    prop_node = docNode.createElement(controlProps{q}); %create entry
    prop_node = docRootNode.appendChild(prop_node); %add entry
    t_form=repmat('%f, ',1,size(controlPropVals{q},2)); t_form=t_form(1:end-2);
    prop_node.appendChild(docNode.createTextNode(sprintf(t_form,controlPropVals{q}))); %append data text child
end

%Adding analysis type
prop_node = docNode.createElement('analysis'); %create entry
prop_node = docRootNode.appendChild(prop_node); %add entry
attr = docNode.createAttribute('type'); %Create id attribute
attr.setNodeValue(FEB_struct.Control.AnalysisType); %Set id text
prop_node.setAttributeNode(attr); %Add id attribute

%Time stepper properties
prop_node = docNode.createElement('time_stepper'); %create entry
prop_node = docRootNode.appendChild(prop_node); %add entry

controlProps=FEB_struct.Control.TimeStepperProperties;
controlPropVals=FEB_struct.Control.TimeStepperValues;
for q=1:1:numel(controlProps)
    prop_prop_node = docNode.createElement(controlProps{q}); %create entry
    prop_prop_node = prop_node.appendChild(prop_prop_node); %add entry
    t_form=repmat('%f, ',1,size(controlPropVals{q},2)); t_form=t_form(1:end-2);
    prop_prop_node.appendChild(docNode.createTextNode(sprintf(t_form,controlPropVals{q}))); %append data text child
end

end