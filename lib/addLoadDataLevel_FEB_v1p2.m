function [docNode]=addLoadDataLevel_FEB_v1p2(docNode,FEB_struct)


%% Add load curves

if  isfield(FEB_struct,'LoadData')
    disp('Adding LoadData level');
    febio_spec = docNode.getDocumentElement;
    
    %Check for LoadData level
    if docNode.getElementsByTagName('<LoadData>').getLength==0; %level does not exist yet
        ElementAddNode = docNode.createElement('LoadData');
        febio_spec.appendChild(ElementAddNode);
        LoadDataNode = docNode.getElementsByTagName('LoadData').item(0);
        IsControlNodeNew=1;
    else
        LoadDataNode = docNode.getElementsByTagName('LoadData').item(0);
        IsControlNodeNew=0;
    end
    docRootNode=LoadDataNode;
    
    %Defining load curves
    disp('----> Defining load curves')
    for q=1:1:numel(FEB_struct.LoadData.LoadCurves.id)
        
        loadPoints=FEB_struct.LoadData.LoadCurves.loadPoints{q};
        
        prop_node = docNode.createElement('loadcurve'); %create entry
        prop_node = docRootNode.appendChild(prop_node); %add entry
        
        attr = docNode.createAttribute('id'); %Create attribute
        attr.setNodeValue(num2str(FEB_struct.LoadData.LoadCurves.id(q))); %Set text
        prop_node.setAttributeNode(attr); %Add attribute
        
        attr = docNode.createAttribute('type'); %Create attribute
        attr.setNodeValue(FEB_struct.LoadData.LoadCurves.type{q}); %Set text
        prop_node.setAttributeNode(attr); %Add attribute
        
        for q2=1:1:size(loadPoints,1)
            loadPoint=loadPoints(q2,:);
            prop_prop_node = docNode.createElement('loadpoint'); %create entry
            prop_prop_node = prop_node.appendChild(prop_prop_node); %add entry
            
            t_form=repmat('%f, ',1,size(loadPoint,2)); t_form=t_form(1:end-2);
            prop_prop_node.appendChild(docNode.createTextNode(sprintf(t_form,loadPoint))); %append data text child
        end
    end
    
end
