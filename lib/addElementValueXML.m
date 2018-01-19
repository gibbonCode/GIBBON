function [domNode]=addElementValueXML(domNode,elementNode,elementValue)

if isnumeric(elementValue) %If it is numeric
    if isrounded(elementValue) %If it looks like an integer
        t_form=repmat('%d, ',1,size(elementValue,2)); t_form=t_form(1:end-2);
    else
        t_form=repmat('%6.7e, ',1,size(elementValue,2)); t_form=t_form(1:end-2);
    end
    elementValue=sprintf(t_form,elementValue);
end

if ischar(elementValue)
    elementNode.appendChild(domNode.createTextNode(elementValue)); %append data text child
else
    error('elementValue should be a character string');
end

end