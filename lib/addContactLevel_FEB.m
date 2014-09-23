function [docNode]=addContactLevel_FEB(docNode,FEB_struct)

disp('Adding Contact field')
rootNode = docNode.getDocumentElement;

contactNode = docNode.createElement('Contact');
contactNode = rootNode.appendChild(contactNode);
    

%% Adding contact components
[docNode]=addContactComponents_FEB(docNode,contactNode,FEB_struct);

