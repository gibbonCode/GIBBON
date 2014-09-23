function [docNode]=addGeometryLevel_TXT(docNode,FEB_struct)

%%

disp('Adding Geometry level');

rootNode = docNode.getDocumentElement;
geometryNode = docNode.createElement('Geometry');
geometryNode = rootNode.appendChild(geometryNode);

%%

geometryNode = docNode.getElementsByTagName('Geometry').item(0);
parent_node = docNode.createElement('Nodes');
parent_node = geometryNode.appendChild(parent_node);

xmlwrite(FEB_struct.run_filename,docNode); %Export to text file
[T]=txtfile2cell(FEB_struct.run_filename); %Import back into cell array

T_top=T(1:end-3);

T_end=T(end-1:end);

%NODE FIELD
disp('----> Adding node field')
if FEB_struct.disp_opt==1;
    hw = waitbar(0,'Adding node entries....');
end

T_node=cell(size(FEB_struct.Geometry.Nodes,1)+2,1);
T_node(1,1)={'		<Nodes>'};
n_steps=size(FEB_struct.Geometry.Nodes,1);
for q_n=1:1:n_steps
    T_node(q_n+1,1)={['			<node id="',sprintf('%u',q_n),'">',sprintf('%10.7e, %10.7e, %10.7e',FEB_struct.Geometry.Nodes(q_n,:)),'</node>']};
    if FEB_struct.disp_opt==1 && rem(q_n,round(n_steps/10))==0;
        waitbar(q_n/n_steps);
    end
end
T_node(end,1)={'		</Nodes>'};

if FEB_struct.disp_opt==1;
    close(hw); drawnow; 
end

% ELEMENT FIELD
disp('----> Adding element field')

numElementsPerSet=cellfun(@(x) size(x,1),FEB_struct.Geometry.Elements);
numElements=sum(numElementsPerSet);
startInds=[0 cumsum(numElementsPerSet(1:end-1))]+1;
endInds=cumsum(numElementsPerSet);

T_elem=cell(numElements+2,1);
T_elem(1,1)={'		<Elements>'};

for q_e1=1:1:numel(FEB_struct.Geometry.Elements)
    elementIndices=startInds(q_e1):endInds(q_e1);
    
    E=FEB_struct.Geometry.Elements{q_e1};
    t_form=repmat('   %u,',1,size(E,2)); t_form=t_form(1:end-1); %text format for sprintf for current element type
    E_type=FEB_struct.Geometry.ElementType{q_e1};
    M=FEB_struct.Geometry.ElementMat{q_e1};
    
    if FEB_struct.disp_opt==1;
        hw = waitbar(0,['Adding ',E_type,' element entries....']);
    end
    disp(['----> Adding ',E_type,' element entries....']);
    
    n_steps=size(E,1);
    for q_e2=1:1:n_steps
        e_ind=elementIndices(q_e2);
        
        T_elem(e_ind+1,1)={['			<',E_type,' id="',sprintf('%u',e_ind),'" mat="',sprintf('%u',M(q_e2)),'">',sprintf(t_form,E(q_e2,:)),'</',E_type,'>']};
        
        if FEB_struct.disp_opt==1 && rem(q_e2,round(n_steps/10))==0;
            waitbar(q_e2/n_steps);
        end
    end
    
    if FEB_struct.disp_opt==1;
        close(hw); drawnow; 
    end
end
T_elem(end,1)={'		</Elements>'};

% ELEMENT DATA FIELD
if isfield(FEB_struct.Geometry,'ElementData');
    disp('----> Adding element data field');
    
    numElementsPerSet=cellfun(@(x) size(x,1),FEB_struct.Geometry.Elements);
    numElementsPerSetCumsum=cumsum(numElementsPerSet);
    elementSize2PerSet=cellfun(@(x) size(x,2),FEB_struct.Geometry.Elements);
    
    elementSize2=zeros(sum(numElementsPerSet),1);
    c=1; 
    for q=1:1:numel(FEB_struct.Geometry.Elements) %element sets        
        elementSize2(c:numElementsPerSetCumsum(q))=elementSize2PerSet(q);
        c=numElementsPerSetCumsum(q)+1;
    end
    
    
    if isfield(FEB_struct.Geometry.ElementData,'Thickness');        
        elementIndicesThickness=FEB_struct.Geometry.ElementData.IndicesForThickness;
        disp('----> Thickness data entries found');
    else
        elementIndicesThickness=[];
    end
    
    if isfield(FEB_struct.Geometry.ElementData,'MatAxis')
        elementIndicesMatAxis=FEB_struct.Geometry.ElementData.MatAxis.ElementIndices;
        disp('----> MatAxis data entries found');
    else
        elementIndicesMatAxis=[];
    end
    
    if isfield(FEB_struct.Geometry.ElementData,'Fiber')
        elementIndicesFiber=FEB_struct.Geometry.ElementData.Fiber.ElementIndices;
        disp('----> MatAxis data entries found');
    else
        elementIndicesFiber=[];
    end
    
    elementIndices=unique([elementIndicesThickness(:);elementIndicesMatAxis(:);elementIndicesFiber(:)]);
    elementIndices=elementIndices(:)';

    T_elementDataMiddle=cell((numel(elementIndices)*2+numel(elementIndicesThickness)+4*numel(elementIndicesMatAxis)+numel(elementIndicesFiber)),1);
        
    n_steps=numel(elementIndices);
    if FEB_struct.disp_opt==1;
        hw = waitbar(0,['Adding ElementData entries....']);
    end
    c=1;
    cw=1;
    for q_e=elementIndices
              
        T_elementDataMiddle(c,1)={['         <element id="',sprintf('%u',q_e),'">']};
         
        %Thickness entries
        if ismember(q_e,elementIndicesThickness)   
 
            t_form=repmat(' %6.7e,',1,elementSize2(q_e)); t_form=t_form(1:end-1); %text format for sprintf for current element type
            
            thicknessData=FEB_struct.Geometry.ElementData.Thickness(FEB_struct.Geometry.ElementData.IndicesForThickness==q_e)*ones(1,elementSize2(q_e));
            thicknessEntryText=sprintf(t_form,thicknessData);
             
            c=c+1;
            T_elementDataMiddle(c,1)={['            <thickness> ',thicknessEntryText,'</thickness>']};

        end %IF Thickness
        
        %Mat_axis and fibre entries
        if ismember(q_e,elementIndicesMatAxis)    
            c=c+1;
            T_elementDataMiddle(c,1)={['            <mat_axis>']};
            
            adData=FEB_struct.Geometry.ElementData.MatAxis.Basis(FEB_struct.Geometry.ElementData.MatAxis.ElementIndices==q_e,:,:);
            
            c=c+1;
            aData=adData(:,:,1);
            aText=sprintf('%6.7e, %6.7e, %6.7e',aData);
            T_elementDataMiddle(c,1)={['               <a>',aText,'</a>']};
            
            c=c+1;
            dData=adData(:,:,2);
            dText=sprintf('%6.7e, %6.7e, %6.7e',dData);
            T_elementDataMiddle(c,1)={['               <d>',dText,'</d>']};
            
            c=c+1;
            T_elementDataMiddle(c,1)={['            </mat_axis>']};
            
        end %IF MatAxis
        
        c=c+1;
        T_elementDataMiddle(c,1)={['         </element>']};
        
        if FEB_struct.disp_opt==1 && rem(cw,round(n_steps/10))==0;
            waitbar(cw/n_steps);
        end
        c=c+1; 
        cw=cw+1;
    end %FOR q_e1
    
    if FEB_struct.disp_opt==1;
        close(hw); drawnow;
    end
    
    %Compose elementData cell
    T_elementData(1,1)={'		<ElementData>'};
    T_elementData(end+1:end+numel(T_elementDataMiddle),1)=T_elementDataMiddle;
    T_elementData(end+1,1)={'		</ElementData>'};
else
    T_elementData={};
end %IF ElementData



%Compose text cell
totalTextCell=T_top;
totalTextCell(end+1:end+numel(T_node))=T_node;
totalTextCell(end+1:end+numel(T_elem))=T_elem;
if ~isempty(T_elementData)
    totalTextCell(end+1:end+numel(T_elementData))=T_elementData;
end
totalTextCell(end+1:end+numel(T_end))=T_end;

cell2txtfile(FEB_struct.run_filename,totalTextCell,1); %Export to text file
docNode = xmlread(FEB_struct.run_filename); %Reimport docNode

end