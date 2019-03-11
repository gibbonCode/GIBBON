function [domNode]=addLoadDataLevel_FEB(domNode,FEB_struct)

%% Add load curves

disp('Adding LoadData level');
rootNode = domNode.getDocumentElement;

%Create LoadData level
loadDataNode = domNode.createElement('LoadData');
loadDataNode = rootNode.appendChild(loadDataNode);

%Defining load curves
disp('----> Defining load curves')
for q=1:1:numel(FEB_struct.LoadData.LoadCurves.id)
    
    loadPoints=FEB_struct.LoadData.LoadCurves.loadPoints{q};
    
    prop_node = domNode.createElement('loadcurve'); %create entry
    prop_node = loadDataNode.appendChild(prop_node); %add entry
    
    attr = domNode.createAttribute('id'); %Create attribute
    attr.setNodeValue(num2str(FEB_struct.LoadData.LoadCurves.id(q))); %Set text
    prop_node.setAttributeNode(attr); %Add attribute
    
    if isfield(FEB_struct.LoadData.LoadCurves,'type');
        attr = domNode.createAttribute('type'); %Create attribute
        attr.setNodeValue(FEB_struct.LoadData.LoadCurves.type{q}); %Set text
        prop_node.setAttributeNode(attr); %Add attribute
    end
    
    if isfield(FEB_struct.LoadData.LoadCurves,'extend');
        attr = domNode.createAttribute('extend'); %Create attribute
        attr.setNodeValue(FEB_struct.LoadData.LoadCurves.extend{q}); %Set text
        prop_node.setAttributeNode(attr); %Add attribute
    end
    
    for q2=1:1:size(loadPoints,1)
        loadPoint=loadPoints(q2,:);
        prop_prop_node = domNode.createElement('loadpoint'); %create entry
        prop_prop_node = prop_node.appendChild(prop_prop_node); %add entry
        
        t_form=repmat('%f, ',1,size(loadPoint,2)); t_form=t_form(1:end-2);
        prop_prop_node.appendChild(domNode.createTextNode(sprintf(t_form,loadPoint))); %append data text child
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
% Copyright (C) 2019  Kevin Mattheus Moerman
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
