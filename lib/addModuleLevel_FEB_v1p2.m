function docNode=addModuleLevel_FEB_v1p2(docNode,FEB_struct)

disp('Adding Module level');

if ~isfield(FEB_struct,'Module')
    %Add defaults
    FEB_struct.Module.Type='solid';
end

febio_spec = docNode.getDocumentElement;

ElementAddNode = docNode.createElement('Module');
febio_spec.appendChild(ElementAddNode);
parent_node = docNode.getElementsByTagName('Module').item(0);

attr = docNode.createAttribute('type'); %Create attribute
attr.setNodeValue(FEB_struct.Module.Type); %Set text
parent_node.setAttributeNode(attr); %Add attribute