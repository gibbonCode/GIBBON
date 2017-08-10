function [docNode]=addControlComponents_FEB(docNode,ControlNode,FEB_struct)

docRootNode=ControlNode;

%Setting Control parameters
if isfield(FEB_struct.Control,'Properties')
    controlProps=FEB_struct.Control.Properties;
    controlPropVals=FEB_struct.Control.Values; %default values
    for q=1:1:numel(controlProps)
        prop_node = docNode.createElement(controlProps{q}); %create entry
        prop_node = docRootNode.appendChild(prop_node); %add entry
        t_form=repmat('%f, ',1,size(controlPropVals{q},2)); t_form=t_form(1:end-2);
        prop_node.appendChild(docNode.createTextNode(sprintf(t_form,controlPropVals{q}))); %append data text child
    end
end

%Adding analysis type
if isfield(FEB_struct.Control,'AnalysisType')
    prop_node = docNode.createElement('analysis'); %create entry
    prop_node = docRootNode.appendChild(prop_node); %add entry
    attr = docNode.createAttribute('type'); %Create id attribute
    attr.setNodeValue(FEB_struct.Control.AnalysisType); %Set id text
    prop_node.setAttributeNode(attr); %Add id attribute
end

%Time stepper properties
if isfield(FEB_struct.Control,'TimeStepperProperties')
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

end
 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
