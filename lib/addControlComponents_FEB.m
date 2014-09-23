function [docNode]=addControlComponents_FEB(docNode,ControlNode,FEB_struct)

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
