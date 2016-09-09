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

