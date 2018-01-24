function [domNode]=addElementValueXML(domNode,elementNode,elementValue)

[elementValue]=vec2strIntDouble(elementValue,'%6.7e');

if ischar(elementValue)
    elementNode.appendChild(domNode.createTextNode(elementValue)); %append data text child
else
    error('elementValue should be a character string');
end

end