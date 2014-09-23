function docNode=addControlLevel_FEB_v1p2(docNode,FEB_struct)

febio_spec = docNode.getDocumentElement;

%Check for Control level
if docNode.getElementsByTagName('Control').getLength==0; %level does not exist yet
    ElementAddNode = docNode.createElement('Control');
    febio_spec.appendChild(ElementAddNode);
    ControlNode = docNode.getElementsByTagName('Control').item(0);
    IsControlNodeNew=1;
else
    ControlNode = docNode.getElementsByTagName('Control').item(0);
    IsControlNodeNew=0;
end

%% Adding control components

[docNode]=addControlComponents_FEB_v1p2(docNode,ControlNode,FEB_struct);


