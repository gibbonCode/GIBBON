function [docNode]=addLoadDataLevel_FEB(docNode,FEB_struct)

%% Add load curves

disp('Adding LoadData level');
rootNode = docNode.getDocumentElement;

%Create LoadData level
loadDataNode = docNode.createElement('LoadData');
loadDataNode = rootNode.appendChild(loadDataNode);

%Defining load curves
disp('----> Defining load curves')
for q=1:1:numel(FEB_struct.LoadData.LoadCurves.id)
    
    loadPoints=FEB_struct.LoadData.LoadCurves.loadPoints{q};
    
    prop_node = docNode.createElement('loadcurve'); %create entry
    prop_node = loadDataNode.appendChild(prop_node); %add entry
    
    attr = docNode.createAttribute('id'); %Create attribute
    attr.setNodeValue(num2str(FEB_struct.LoadData.LoadCurves.id(q))); %Set text
    prop_node.setAttributeNode(attr); %Add attribute
    
    if isfield(FEB_struct.LoadData.LoadCurves,'type');
        attr = docNode.createAttribute('type'); %Create attribute
        attr.setNodeValue(FEB_struct.LoadData.LoadCurves.type{q}); %Set text
        prop_node.setAttributeNode(attr); %Add attribute
    end
    
    if isfield(FEB_struct.LoadData.LoadCurves,'extend');
        attr = docNode.createAttribute('extend'); %Create attribute
        attr.setNodeValue(FEB_struct.LoadData.LoadCurves.extend{q}); %Set text
        prop_node.setAttributeNode(attr); %Add attribute
    end
    
    for q2=1:1:size(loadPoints,1)
        loadPoint=loadPoints(q2,:);
        prop_prop_node = docNode.createElement('loadpoint'); %create entry
        prop_prop_node = prop_node.appendChild(prop_prop_node); %add entry
        
        t_form=repmat('%f, ',1,size(loadPoint,2)); t_form=t_form(1:end-2);
        prop_prop_node.appendChild(docNode.createTextNode(sprintf(t_form,loadPoint))); %append data text child
    end
end
 
%% 
% ********** _license boilerplate_ **********
% 
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
