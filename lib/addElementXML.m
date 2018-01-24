function [varargout]=addElementXML(domNode,rootNode,elementName,elementValue)

try
    elementNode = domNode.createElement(elementName); %create entry
    elementNode = rootNode.appendChild(elementNode); %add entry
catch
    % Get correct rootNode
    nodeList=domNode.getElementsByTagName(rootNode.getNodeName);
    rootNode=nodeList.item(nodeList.getLength-1);
    
    elementNode = domNode.createElement(elementName); %create entry
    elementNode = rootNode.appendChild(elementNode); %add entry    
end

if ~isempty(elementValue)
    [elementValue]=vec2strIntDouble(elementValue,'%6.7e');
    elementNode.appendChild(domNode.createTextNode(elementValue)); %append data text child
end

%Collect output
varargout{1}=domNode;
varargout{2}=elementNode;
varargout{3}=rootNode;

end