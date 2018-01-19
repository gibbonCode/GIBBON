function [domNode]=addAttributeXML(domNode,elementNode,attrName,attrValue)

[attrValue]=vec2strIntDouble(attrValue,'%6.7e');

if ischar(attrValue)
    attr = domNode.createAttribute(attrName); %Create attribute
    attr.setNodeValue(attrValue); %Set text
    elementNode.setAttributeNode(attr); %Add attribute
else
    error('attrValue should be a character string');
end

end