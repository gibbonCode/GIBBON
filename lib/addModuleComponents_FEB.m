function [docNode]=addModuleComponents_FEB(docNode,moduleNode,FEB_struct)

%Check for parameters or use default
if ~isfield(FEB_struct,'Module')
    %Add defaults
    FEB_struct.Module.Type='solid';
end

attr = docNode.createAttribute('type'); %Create attribute
attr.setNodeValue(FEB_struct.Module.Type); %Set text
moduleNode.setAttributeNode(attr); %Add attribute

