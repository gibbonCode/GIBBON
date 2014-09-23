function docNode=addControlLevel_FEB(docNode,FEB_struct)

disp('Adding Control level')

%% Create control level
rootNode = docNode.getDocumentElement;
controlNode = docNode.createElement('Control');
controlNode = rootNode.appendChild(controlNode);

%% Adding control components
[docNode]=addControlComponents_FEB(docNode,controlNode,FEB_struct);

end
