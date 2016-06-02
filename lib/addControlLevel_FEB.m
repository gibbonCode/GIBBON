function domNode=addControlLevel_FEB(domNode,FEB_struct)

% function [domNode]=addControlLevel_FEB(domNode,FEB_struct)
% ------------------------------------------------------------------------
% Adds control information to the XML object domNode based on
% the FEBio structure FEB_struct. 
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2016/06/02
%------------------------------------------------------------------------

%%

disp('Adding Control level')

%% Create control level
rootNode = domNode.getDocumentElement;
controlNode = domNode.createElement('Control');
controlNode = rootNode.appendChild(controlNode);

%% Adding control components
[domNode]=addControlComponents_FEB(domNode,controlNode,FEB_struct);

end
