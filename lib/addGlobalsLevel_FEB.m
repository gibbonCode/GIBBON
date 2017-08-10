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
