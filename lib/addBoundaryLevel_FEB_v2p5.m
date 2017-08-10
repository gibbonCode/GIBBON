function [domNode]=addBoundaryLevel_FEB_v2p5(domNode,FEB_struct)

% function [domNode]=addBoundaryLevel_FEB_v2p5(domNode,FEB_struct)
% ------------------------------------------------------------------------
% Adds boundary condition information to the XML object domNode based on
% the FEBio structure FEB_struct. 
% 
% Fixed and prescribed boundary conditions are supported
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2016/06/02
%------------------------------------------------------------------------

%%

disp('Adding Boundary level')

rootNode = domNode.getDocumentElement;

boundaryNode = domNode.createElement('Boundary');
boundaryNode = rootNode.appendChild(boundaryNode);

%% Adding boundary components

[domNode]=addBoundaryComponents_FEB(domNode,boundaryNode,FEB_struct);

 
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
