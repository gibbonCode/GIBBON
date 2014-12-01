function docNode=addModuleLevel_FEB(docNode,FEB_struct)

disp('Adding Module level');

%% Create Module level
rootNode = docNode.getDocumentElement;
moduleNode = docNode.createElement('Module');
moduleNode = rootNode.appendChild(moduleNode);

%% Adding boundary components

[docNode]=addModuleComponents_FEB(docNode,moduleNode,FEB_struct);


