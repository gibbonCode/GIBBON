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
% 
%     GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
%     image segmentation, image-based modeling, meshing, and finite element
%     analysis.
% 
%     Copyright (C) 2017  Kevin Mattheus Moerman
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
