function [domNode]=addElementSetXML(domNode,rootNode,parseStruct)

fieldNameSet = fieldnames(parseStruct);
for q_field=1:1:numel(fieldNameSet)
    currentElementName=fieldNameSet{q_field};
    currentElementValue=parseStruct.(fieldNameSet{q_field});    
    [domNode]=addElementXML(domNode,rootNode,currentElementName,currentElementValue);    
end

end