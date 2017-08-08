function [domNode]=addContactLevel_FEB(domNode,FEB_struct)

% function [domNode]=addContactLevel_FEB(domNode,FEB_struct)
% ------------------------------------------------------------------------
% Adds contact information to the XML object domNode based on
% the FEBio structure FEB_struct. 
% 
% Sticky and sliding contact formulations are supported
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2016/06/02
%------------------------------------------------------------------------

disp('Adding Contact field')
rootNode = domNode.getDocumentElement;

contactNode = domNode.createElement('Contact');
contactNode = rootNode.appendChild(contactNode);
    

%% Adding contact components
[domNode]=addContactComponents_FEB(domNode,contactNode,FEB_struct);

 
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
