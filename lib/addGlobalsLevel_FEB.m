function domNode=addGlobalsLevel_FEB(domNode,FEB_struct)

% function [domNode]=addGlobalsLevel_FEB(domNode,FEB_struct)
% ------------------------------------------------------------------------
% Adds globals section to the XML object domNode based on
% the FEBio structure FEB_struct. 
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2016/06/02
%------------------------------------------------------------------------
%%

disp('Adding Globals level');

%% 
if ~isfield(FEB_struct,'Globals')
    %Add defaults
    FEB_struct.Globals.Constants.Names={'T','R','Fc'};
    FEB_struct.Globals.Constants.Entries={0,0,0};
end

rootNode = domNode.getDocumentElement;

%% Create Globals level

globalsNode = domNode.createElement('Globals');
globalsNode = rootNode.appendChild(globalsNode);

% Adding Constants
if isfield(FEB_struct.Globals,'Constants')    
    %Adding constants node
    constants_node = domNode.createElement('Constants'); %create material entry
    constants_node = globalsNode.appendChild(constants_node); %add material entry
    for q=1:1:numel(FEB_struct.Globals.Constants.Names)
        constantName=FEB_struct.Globals.Constants.Names{q}; %Constant name
        constantEntry=FEB_struct.Globals.Constants.Entries{q}; %Constant entry      
        attribute_node = domNode.createElement(constantName); %create entry
        attribute_node = constants_node.appendChild(attribute_node); %add entry
        if ischar(constantEntry)
            attribute_node.appendChild(domNode.createTextNode(constantEntry)); %append data text child
        else
            t_form=repmat('%6.7e, ',1,size(constantEntry,2)); t_form=t_form(1:end-2);
            attribute_node.appendChild(domNode.createTextNode(sprintf(t_form,constantEntry))); %append data text child
        end        
    end    
end

% % Adding Generations
% if isfield(FEB_struct.Globals,'Generations')    
%     %Adding constants node
%     generations_node = docNode.createElement('Generations'); %create material entry
%     generations_node = globalsNode.appendChild(generations_node); %add material entry
%     for q=1:1:numel(FEB_struct.Globals.Generations.id)
%         attribute_node = docNode.createElement('gen'); %create entry
%         attribute_node = generations_node.appendChild(attribute_node); %add entry
%         
%         attr = docNode.createAttribute('id'); %Create attribute
%         attr.setNodeValue(num2str(FEB_struct.Globals.Generations.id(q))); %Set text
%         attribute_node.setAttributeNode(attr); %Add attribute
%         
%         tGammaEntry=FEB_struct.Globals.Generations.tGamma(q);
%         t_form=repmat('%6.7e, ',1,size(tGammaEntry,2)); t_form=t_form(1:end-2);
%         attribute_node.appendChild(docNode.createTextNode(sprintf(t_form,tGammaEntry))); %append data text child        
%     end    
% end

 
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
