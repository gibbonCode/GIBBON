function docNode=addGeometryLevel_FEB(docNode,FEB_struct)

%%

disp('Adding Geometry level');

rootNode = docNode.getDocumentElement;
geometryNode = docNode.createElement('Geometry');
geometryNode = rootNode.appendChild(geometryNode);

%% Adding node field

disp('----> Adding node field');

parent_node = docNode.createElement('Nodes');
parent_node = geometryNode.appendChild(parent_node);

if FEB_struct.disp_opt==1;
    hw = waitbar(0,'Adding node entries....');
end
n_steps=size(FEB_struct.Geometry.Nodes,1);
for q_n=1:1:n_steps
    node_node = docNode.createElement('node'); %create node entry
    node_node = parent_node.appendChild(node_node); %add node entry
    attr = docNode.createAttribute('id'); %Create id attribute
    attr.setNodeValue(sprintf('%u',q_n)); %Set id text
    node_node.setAttributeNode(attr); %Add id attribute
    node_node.appendChild(docNode.createTextNode(sprintf('%6.7e, %6.7e, %6.7e',FEB_struct.Geometry.Nodes(q_n,:)))); %append data text child
    if FEB_struct.disp_opt==1 && rem(q_n,round(n_steps/10))==0;
        waitbar(q_n/n_steps);
    end
end

if FEB_struct.disp_opt==1;
    close(hw); drawnow;
end

%% Adding element field

disp('----> Adding element field')

geometryNode = docNode.getElementsByTagName('Geometry').item(0); %Get GEONode

numElementsPerSet=cellfun(@(x) size(x,1),FEB_struct.Geometry.Elements);
startInds=[0 cumsum(numElementsPerSet(1:end-1))]+1;
endInds=cumsum(numElementsPerSet);

for q_e1=1:1:numel(FEB_struct.Geometry.Elements) %for all element sets
    
    E_type=FEB_struct.Geometry.ElementType{q_e1}; %Element type for current set
    partName=FEB_struct.Geometry.ElementsPartName{q_e1}; %Element set part name
    
    if FEB_struct.disp_opt==1;
        hw = waitbar(0,['Adding ',E_type,' element entries....']);
    end
    disp(['----> Adding ',E_type,' element entries....']);
    
    elementIndices=startInds(q_e1):endInds(q_e1);
    
    E=FEB_struct.Geometry.Elements{q_e1}; %The current element set
    
    t_form=repmat('   %u,',1,size(E,2)); t_form=t_form(1:end-1); %text format for sprintf for current element type
    
    M=FEB_struct.Geometry.ElementMat{q_e1}; %Current elements set's material indices
    matIdSet=unique(M(:)); %Set of unique material indices
    numMat=numel(matIdSet);
    
    for q_mat=1:1:numMat
        
        matId=matIdSet(q_mat);
        logicMat=M==matId; %Logic for all elements for matId
        E_mat=E(logicMat,:);
        elementIndices_mat=elementIndices(logicMat);
        
        if numMat>1
            partNameSet=[partName,'_',sprintf('%u',matId)];
        else
            partNameSet=partName;
        end
        
        %Add Elements field
        parent_node = docNode.createElement('Elements');
        parent_node = geometryNode.appendChild(parent_node);
        
        %Set type attribute
        attr = docNode.createAttribute('type'); %Create id attribute
        attr.setNodeValue(E_type); %Set id text
        parent_node.setAttributeNode(attr); %Add id attribute
        
        %Set mat attribute
        attr = docNode.createAttribute('mat'); %Create id attribute
        attr.setNodeValue(sprintf('%u',matId)); %Set id text
        parent_node.setAttributeNode(attr); %Add id attribute
        
        %Set mat attribute
        attr = docNode.createAttribute('name'); %Create id attribute
        attr.setNodeValue(partNameSet); %Set id text
        parent_node.setAttributeNode(attr); %Add id attribute
        
        n_steps=size(E_mat,1);
        for q_e2=1:1:n_steps
            e_ind=elementIndices_mat(q_e2);
            
            element_node= docNode.createElement('elem');
            element_node = parent_node.appendChild(element_node);
            
            attr = docNode.createAttribute('id'); %Create id attribute
            attr.setNodeValue(sprintf('%u',e_ind)); %Set id text
            element_node.setAttributeNode(attr); %Add id attribute
            
            element_node.appendChild(docNode.createTextNode(sprintf(t_form,E_mat(q_e2,:)))); %append data text child
            
            if FEB_struct.disp_opt==1 && rem(q_e2,round(n_steps/10))==0;
                waitbar(q_e2/n_steps);
            end
        end
        
    end
    
    if FEB_struct.disp_opt==1;
        close(hw); drawnow;
    end
end

%% ElementData for thickness

if isfield(FEB_struct.Geometry,'ElementData');
    if isfield(FEB_struct.Geometry.ElementData,'Thickness');
        geometryNode = docNode.getElementsByTagName('Geometry').item(0); %Get GEONode
        
        numElementsPerSet=cellfun(@(x) size(x,1),FEB_struct.Geometry.Elements);
        startInds=[0 cumsum(numElementsPerSet(1:end-1))]+1;
        endInds=cumsum(numElementsPerSet);
        
        addedElementDataField=0;
        for q_e1=1:1:numel(FEB_struct.Geometry.Elements)
            elementIndices=startInds(q_e1):endInds(q_e1);
            logicThicknessSet=ismember(elementIndices,FEB_struct.Geometry.ElementData.IndicesForThickness);
            
            if any(logicThicknessSet) %If any thicknesses are to be set for the current element group
                
                E=FEB_struct.Geometry.Elements{q_e1};
                E_type=FEB_struct.Geometry.ElementType{q_e1};
                
                if FEB_struct.disp_opt==1;
                    hw = waitbar(0,['Adding ',E_type,' element thickness data entries....']);
                end
                c=1;
                
                disp(['----> Adding ',E_type,' element thickness data entries....']);
                n_steps=nnz(logicThicknessSet);
                for e_ind=elementIndices(logicThicknessSet)
                    
                    if addedElementDataField==0
                        parent_node1 = docNode.createElement('ElementData');
                        parent_node1 = geometryNode.appendChild(parent_node1);
                        addedElementDataField=1;
                    end
                    
                    %Creating ElementData entries
                    element_node = docNode.createElement('element');
                    element_node = parent_node1.appendChild(element_node);
                    attr = docNode.createAttribute('id');
                    attr.setNodeValue(sprintf('%u',e_ind));
                    element_node.setAttributeNode(attr);
                    
                    %Adding thickness level
                    t_form=repmat(' %6.7e,',1,size(E,2)); t_form=t_form(1:end-1); %text format for sprintf for current element type
                    
                    v_text=sprintf(t_form,FEB_struct.Geometry.ElementData.Thickness(FEB_struct.Geometry.ElementData.IndicesForThickness==e_ind)*ones(1,size(E,2)));
                    thickness_node = docNode.createElement('thickness');
                    element_node.appendChild(thickness_node);
                    thickness_node.appendChild(docNode.createTextNode(v_text)); %Adding thickness data text
                    
                    if FEB_struct.disp_opt==1 && rem(c,round(n_steps/10))==0;
                        waitbar(c/n_steps);
                        
                    end
                    c=c+1;
                end
                if FEB_struct.disp_opt==1;
                    close(hw); drawnow;
                end
            end
        end
    end
end

%% ElementData for fibre directions and MatAxis

% Checking for fiber / mat_axis information
if isfield(FEB_struct.Geometry,'ElementData') %If ElementData exists
    if isfield(FEB_struct.Geometry.ElementData,'MatAxis') %If MatAxis data exists
        docNode=addMatAxisFibreElementData(docNode,FEB_struct);
    end
end

%% Adding surface field
if isfield(FEB_struct.Geometry,'Surface');
    disp('----> Adding surface field');
    
    for q_set=1:1:numel(FEB_struct.Geometry.Surface)
        
        F=FEB_struct.Geometry.Surface{q_set}.Set; %Faces
        surfaceType=FEB_struct.Geometry.Surface{q_set}.Type; %Surface type
        
        
        %Create surface level
        parent_node = docNode.createElement('Surface');
        parent_node = geometryNode.appendChild(parent_node);
        
        if ~isfield(FEB_struct.Geometry.Surface{q_set},'Name')
            FEB_struct.Geometry.Surface{q_set}.Name=['Surface_',sprintf('%u',q_set)]; %Surface type
        end
        surfaceName=FEB_struct.Geometry.Surface{q_set}.Name; %Surface type
        
        attr = docNode.createAttribute('name'); %Create id attribute
        attr.setNodeValue(surfaceName); %Set id text
        parent_node.setAttributeNode(attr); %Add id attribute
        
        numFaces=size(F,1);
        for q_face=1:1:numFaces %For all faces
            
            %Add element entry
            element_node = docNode.createElement(surfaceType);
            element_node = parent_node.appendChild(element_node);
            
            %Add id attribute
            attr = docNode.createAttribute('id'); %Create id attribute
            attr.setNodeValue(sprintf('%u',q_face)); %Set id text
            element_node.setAttributeNode(attr); %Add id attribute
            
            %Add node indices for surface element
            t_form=repmat('   %u,',1,size(F,2)); t_form=t_form(1:end-1); %text format for sprintf for current element type
            element_node.appendChild(docNode.createTextNode(sprintf(t_form,F(q_face,:)))); %append data text child
            
        end
    end
end

%% Adding NodeSet section

% <NodeSet name="nodeset1">
%       <node id="1"/>
%       <node id="2"/>
%       <node id="101"/>
%       <node id="102"/>
% </NodeSet>

writeNodeSetType=1; %Override as 2 if the old format is desired which has the form form: <NodeSet name="bcSupportList_Z">  1,  2,  3,  4,  5,  6 </NodeSet>

if isfield(FEB_struct.Geometry,'NodeSet');
    disp('----> Adding NodeSet field');
    
    for q_set=1:1:numel(FEB_struct.Geometry.NodeSet)
        
        nodeSetIndices=FEB_struct.Geometry.NodeSet{q_set}.Set; %Node indices
        nodeSetIndices=nodeSetIndices(:)'; %Row vector
        
        %Create NodeSet level
        parent_node = docNode.createElement('NodeSet');
        parent_node = geometryNode.appendChild(parent_node);
        
        if ~isfield(FEB_struct.Geometry.NodeSet{q_set},'Name')
            FEB_struct.Geometry.NodeSet{q_set}.Name=['NodeSet_',sprintf('%u',q_set)]; %Surface type
        end
        nodeSetName=FEB_struct.Geometry.NodeSet{q_set}.Name; %Node set name
        
        attr = docNode.createAttribute('name'); %Create id attribute
        attr.setNodeValue(nodeSetName); %Set id text
        parent_node.setAttributeNode(attr); %Add id attribute
        
        switch writeNodeSetType
            case 1 %New since FEBio 2.2
                numNodes=numel(nodeSetIndices);
                for q_node=1:1:numNodes %For all nodes
                    
                    %Add element entry
                    element_node = docNode.createElement('node');
                    element_node = parent_node.appendChild(element_node);
                    
                    %Add id attribute
                    attr = docNode.createAttribute('id'); %Create id attribute
                    attr.setNodeValue(sprintf('%u',q_node)); %Set id text
                    element_node.setAttributeNode(attr); %Add id attribute
                end
            case 2 %OLD : Will be removed in future releases
                %Add node index text field
                maxWidth=25;
                maxCharLength=numel(sprintf('%u',max(nodeSetIndices)));
                t_form_sub=['%',sprintf('%u',maxCharLength+1),'u,'];
                t_form=[repmat(t_form_sub,1,maxWidth),'\n'];
                textSet=sprintf(t_form,nodeSetIndices);
                if mod(numel(nodeSetIndices),maxWidth)==0 %Devisable by 5 so a , and \n appear at the end
                    textSet=textSet(1:end-2);%Take away last \n and ,
                elseif mod(numel(nodeSetIndices),maxWidth)~=0 %A remainder leading to only an extra ,
                    textSet=textSet(1:end-1); %Take away last ,
                end
                parent_node.appendChild(docNode.createTextNode(textSet)); %append data text child
        end
        
    end    
end

end



function docNode=addMatAxisFibreElementData(docNode,FEB_struct)

Vf=FEB_struct.Geometry.ElementData.MatAxis.Basis;
E_ind=FEB_struct.Geometry.ElementData.MatAxis.ElementIndices;

%Get Geometry level
GEONode = docNode.getElementsByTagName('Geometry').item(0);
docRootNode =  GEONode;

%Check for ElementData level
if GEONode.getElementsByTagName('ElementData').getLength==0; %ElementData level does not exist yet
    %Adding ElementData level
    parent_node = docNode.createElement('ElementData');
    parent_node = docRootNode.appendChild(parent_node);
    new_entry=1;
else %ElementData level already exists
    %Finding existing element entries
    no_entries=docRootNode.getElementsByTagName('element').getLength;
    ID_no=zeros(1,no_entries);
    for i=0:1:no_entries-1
        ID_no(i+1)=str2double(docRootNode.getElementsByTagName('element').item(i).getAttribute('id').toCharArray()');
    end
    new_entry=0;
end
docRootNode =  GEONode.getElementsByTagName('ElementData').item(0);

%Creating ElementData entries
if FEB_struct.disp_opt==1;
    hw = waitbar(0,'Creating MatAxis entries...');
end
disp('----> Creating MatAxis entries')

n_steps = numel(E_ind);
for i=1:n_steps
    if new_entry==1 || ~any(ID_no==E_ind(i)); %Need to add element level
        element_node = docNode.createElement('element');
        element_node = docRootNode.appendChild(element_node);
        attr = docNode.createAttribute('id');
        attr.setNodeValue(sprintf('%u',E_ind(i)));
        element_node.setAttributeNode(attr);
    else %element level already exists need to check for fiber entries
        element_node=docRootNode.getElementsByTagName('element').item(find(ID_no==E_ind(i))-1);
    end
    
    if size(Vf,3)==1
        %case 'transiso'
        v_text=sprintf('%6.7e, %6.7e, %6.7e',Vf(i,:));
        if element_node.getElementsByTagName('fiber').getLength==0
            %Adding fiber level
            fiber_node = docNode.createElement('fiber');
            element_node.appendChild(fiber_node);
            %Adding fiber data text
            fiber_node.appendChild(docNode.createTextNode(v_text));
        else %Fiber entries exist, overwriting existing
            element_node.getElementsByTagName('fiber').item(0).getFirstChild.setData(char(v_text));
        end
    else
        %case 'ortho'
        a_text=sprintf('%6.7e, %6.7e, %6.7e',Vf(i,:,1));
        d_text=sprintf('%6.7e, %6.7e, %6.7e',Vf(i,:,2));
        
        if element_node.getElementsByTagName('mat_axis').getLength==0
            %Adding mat_axis level
            mat_axis_node = docNode.createElement('mat_axis');
            element_node.appendChild(mat_axis_node);
            %Adding a level
            a_node = docNode.createElement('a');
            mat_axis_node.appendChild(a_node);
            %Adding a data text
            a_node.appendChild(docNode.createTextNode(a_text));
            %Adding d level
            d_node = docNode.createElement('d');
            mat_axis_node.appendChild(d_node);
            %Adding d data text
            d_node.appendChild(docNode.createTextNode(d_text));
        else %Mat_axis entries exist, overwriting existing
            %Finding mat_axis level
            mat_axis_node=element_node.getElementsByTagName('mat_axis').item(0);
            %Overwriting a level
            mat_axis_node.getElementsByTagName('a').item(0).getFirstChild.setData(char(a_text));
            %Overwriting d level
            mat_axis_node.getElementsByTagName('d').item(0).getFirstChild.setData(char(d_text));
        end
    end
    if FEB_struct.disp_opt==1 && rem(i,round(n_steps/10))==0;
        waitbar(i/n_steps);
    end
end

if FEB_struct.disp_opt==1;
    close(hw); drawnow;
end

end
