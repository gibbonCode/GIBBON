function [domNode]=addAttributeSetXML(domNode,elementNode,parseStruct)

fieldNameSet = fieldnames(parseStruct);
for q_field=1:1:numel(fieldNameSet)
    attrName=fieldNameSet{q_field};
    attrValue=parseStruct.(attrName);
    [domNode]=addAttributeXML(domNode,elementNode,attrName,attrValue);
end

