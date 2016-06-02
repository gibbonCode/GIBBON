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

