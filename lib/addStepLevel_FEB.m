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
