function [domNode]=addStepLevel_FEB(domNode,FEB_struct)

rootNode = domNode.getDocumentElement;

for q_step=1:1:numel(FEB_struct.Step)
    
    % Adding step level
    disp('Adding Step level');
    stepNode = domNode.createElement('Step');
    rootNode.appendChild(stepNode);
    
    attr = domNode.createAttribute('name'); %Create id attribute
    attr.setNodeValue(['Step_',sprintf('%u',q_step)]); %Set text
    stepNode.setAttributeNode(attr); %Add attribute
        
    % Adding module components
    disp('----> Adding Module field');
    moduleNode = domNode.createElement('Module'); %create entry
    moduleNode = stepNode.appendChild(moduleNode); %add entry
    FEB_struct_temp=FEB_struct.Step{q_step};
    FEB_struct_temp.disp_opt=FEB_struct.disp_opt;
    [domNode]=addModuleComponents_FEB(domNode,moduleNode,FEB_struct_temp);
    
    % Adding boundary components
    if isfield(FEB_struct.Step{q_step},'Boundary')
        disp('----> Adding Boundary field');        
        boundaryNode = domNode.createElement('Boundary'); %create entry
        boundaryNode = stepNode.appendChild(boundaryNode); %add entry
        FEB_struct_temp=FEB_struct.Step{q_step};
        FEB_struct_temp.disp_opt=FEB_struct.disp_opt;
        [domNode]=addBoundaryComponents_FEB(domNode,boundaryNode,FEB_struct_temp);
    end
    
    %Adding contact components
    if isfield(FEB_struct.Step{q_step},'Contact')
        disp('----> Adding Contact field');
        contactNode = domNode.createElement('Contact'); %create entry
        contactNode = stepNode.appendChild(contactNode); %add entry
        FEB_struct_temp=FEB_struct.Step{q_step};
        FEB_struct_temp.disp_opt=FEB_struct.disp_opt;
        [domNode]=addContactComponents_FEB(domNode,contactNode,FEB_struct_temp);
    end
    
    % Adding control components
    if isfield(FEB_struct.Step{q_step},'Control')
        disp('----> Adding Control field');
        controlNode = domNode.createElement('Control'); %create entry
        controlNode = stepNode.appendChild(controlNode); %add entry
        FEB_struct_temp=FEB_struct.Step{q_step};
        FEB_struct_temp.disp_opt=FEB_struct.disp_opt;
        [domNode]=addControlComponents_FEB(domNode,controlNode,FEB_struct_temp);
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
% Copyright (C) 2017  Kevin Mattheus Moerman
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
