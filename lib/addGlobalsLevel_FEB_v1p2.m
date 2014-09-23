function docNode=addGlobalsLevel_FEB_v1p2(docNode,FEB_struct)

disp('Adding Globals level');

if ~isfield(FEB_struct,'Globals')
    %Add defaults
    FEB_struct.Globals.Constants.Names={'T','R','Fc'};
    FEB_struct.Globals.Constants.Entries={0,0,0};
end

febio_spec = docNode.getDocumentElement;

ElementAddNode = docNode.createElement('Globals');
febio_spec.appendChild(ElementAddNode);
parent_node = docNode.getElementsByTagName('Globals').item(0);

% Adding Constants
if isfield(FEB_struct.Globals,'Constants')    
    %Adding constants node
    constants_node = docNode.createElement('Constants'); %create material entry
    constants_node = parent_node.appendChild(constants_node); %add material entry
    for q=1:1:numel(FEB_struct.Globals.Constants.Names)
        constantName=FEB_struct.Globals.Constants.Names{q}; %Constant name
        constantEntry=FEB_struct.Globals.Constants.Entries{q}; %Constant entry      
        attribute_node = docNode.createElement(constantName); %create entry
        attribute_node = constants_node.appendChild(attribute_node); %add entry
        if ischar(constantEntry)
            attribute_node.appendChild(docNode.createTextNode(constantEntry)); %append data text child
        else
            t_form=repmat('%6.7e, ',1,size(constantEntry,2)); t_form=t_form(1:end-2);
            attribute_node.appendChild(docNode.createTextNode(sprintf(t_form,constantEntry))); %append data text child
        end        
    end    
end

% Adding Generations
if isfield(FEB_struct.Globals,'Generations')    
    %Adding constants node
    generations_node = docNode.createElement('Generations'); %create material entry
    generations_node = parent_node.appendChild(generations_node); %add material entry
    for q=1:1:numel(FEB_struct.Globals.Generations.id)
        attribute_node = docNode.createElement('gen'); %create entry
        attribute_node = generations_node.appendChild(attribute_node); %add entry
        
        attr = docNode.createAttribute('id'); %Create attribute
        attr.setNodeValue(num2str(FEB_struct.Globals.Generations.id(q))); %Set text
        attribute_node.setAttributeNode(attr); %Add attribute
        
        tGammaEntry=FEB_struct.Globals.Generations.tGamma(q);
        t_form=repmat('%6.7e, ',1,size(tGammaEntry,2)); t_form=t_form(1:end-2);
        attribute_node.appendChild(docNode.createTextNode(sprintf(t_form,tGammaEntry))); %append data text child        
    end    
end
