function [docNode]=addStepLevel_FEB(docNode,FEB_struct)

rootNode = docNode.getDocumentElement;

for q_step=1:1:numel(FEB_struct.Step)
    
    % Adding step level
    disp('Adding Step level');
    stepNode = docNode.createElement('Step');
    rootNode.appendChild(stepNode);
    
    attr = docNode.createAttribute('name'); %Create id attribute
    attr.setNodeValue(['Step_',sprintf('%u',q_step)]); %Set text
    stepNode.setAttributeNode(attr); %Add attribute
        
    % Adding module components
    disp('----> Adding Module field');
    moduleNode = docNode.createElement('Module'); %create entry
    moduleNode = stepNode.appendChild(moduleNode); %add entry
    FEB_struct_temp=FEB_struct.Step{q_step};
    FEB_struct_temp.disp_opt=FEB_struct.disp_opt;
    [docNode]=addModuleComponents_FEB(docNode,moduleNode,FEB_struct_temp);
    
    % Adding boundary components
    if isfield(FEB_struct.Step{q_step},'Boundary')
        disp('----> Adding Boundary field');        
        boundaryNode = docNode.createElement('Boundary'); %create entry
        boundaryNode = stepNode.appendChild(boundaryNode); %add entry
        FEB_struct_temp=FEB_struct.Step{q_step};
        FEB_struct_temp.disp_opt=FEB_struct.disp_opt;
        [docNode]=addBoundaryComponents_FEB(docNode,boundaryNode,FEB_struct_temp);
    end
    
    %Adding contact components
    if isfield(FEB_struct.Step{q_step},'Contact')
        disp('----> Adding Contact field');
        contactNode = docNode.createElement('Contact'); %create entry
        contactNode = stepNode.appendChild(contactNode); %add entry
        FEB_struct_temp=FEB_struct.Step{q_step};
        FEB_struct_temp.disp_opt=FEB_struct.disp_opt;
        [docNode]=addContactComponents_FEB(docNode,contactNode,FEB_struct_temp);
    end
    
    % Adding control components
    if isfield(FEB_struct.Step{q_step},'Control')
        disp('----> Adding Control field');
        controlNode = docNode.createElement('Control'); %create entry
        controlNode = stepNode.appendChild(controlNode); %add entry
        FEB_struct_temp=FEB_struct.Step{q_step};
        FEB_struct_temp.disp_opt=FEB_struct.disp_opt;
        [docNode]=addControlComponents_FEB(docNode,controlNode,FEB_struct_temp);
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
