function [varargout]=addElementXML(domNode,rootNode,elementName,elementValue)

elementNode = domNode.createElement(elementName); %create entry
elementNode = rootNode.appendChild(elementNode); %add entry

if ~isempty(elementValue)
    [elementValue]=vec2strIntDouble(elementValue,'%6.7e');
    elementNode.appendChild(domNode.createTextNode(elementValue)); %append data text child
end

varargout{1}=domNode;
varargout{2}=elementNode;

end