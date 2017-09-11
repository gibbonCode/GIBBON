function [domNode]=addGeometryLevel_TXT(domNode,FEB_struct)

% function [domNode]=addGeometryLevel_TXT(domNode,FEB_struct)
% ------------------------------------------------------------------------
% Adds geometry information to the XML object domNode based on
% the FEBio structure FEB_struct. This functions avoids most of the MATLAB
% XML processing and focusses in stead on text file operations which are
% faster for very large XML files. 
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2016/06/02
%------------------------------------------------------------------------

%% Set default display setting if missing
if ~isfield(FEB_struct,'disp_opt')
    FEB_struct.disp_opt=0;
end
%%

disp('Adding Geometry level');

rootNode = domNode.getDocumentElement;
geometryNode = domNode.createElement('Geometry');
geometryNode = rootNode.appendChild(geometryNode);

%% Adding node field

disp('----> Adding node field');

parent_node = domNode.createElement('Nodes');
parent_node = geometryNode.appendChild(parent_node);

xmlwrite(FEB_struct.run_filename,domNode); %Export to text file
[T]=txtfile2cell(FEB_struct.run_filename); %Import back into cell array

T_top=T(1:end-3);
T_end=T(end-1:end);

if FEB_struct.disp_opt==1
    hw = waitbar(0,'Adding node entries....');
end

T_node=cell(size(FEB_struct.Geometry.Nodes,1)+2,1);
T_node(1,1)={'		<Nodes>'};
n_steps=size(FEB_struct.Geometry.Nodes,1);
for q_n=1:1:n_steps
    T_node(q_n+1,1)={['			<node id="',sprintf('%u',q_n),'">',sprintf('%10.7e, %10.7e, %10.7e',FEB_struct.Geometry.Nodes(q_n,:)),'</node>']};
    if FEB_struct.disp_opt==1 && rem(q_n,round(n_steps/10))==0
        waitbar(q_n/n_steps);
    end
end
T_node(end,1)={'		</Nodes>'};

if FEB_struct.disp_opt==1
    close(hw); drawnow;
end

%% Adding element field

disp('----> Adding element field')

numElementsPerSet=cellfun(@(x) size(x,1),FEB_struct.Geometry.Elements);
startInds=[0 cumsum(numElementsPerSet(1:end-1))]+1;
endInds=cumsum(numElementsPerSet);
T_elem={};
for q_e1=1:1:numel(FEB_struct.Geometry.Elements) %for all element sets
    
    E_type=FEB_struct.Geometry.ElementType{q_e1}; %Element type for current set
    partName=FEB_struct.Geometry.ElementsPartName{q_e1}; %Element set part name
    
    if FEB_struct.disp_opt==1
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
        
        n_steps=size(E_mat,1);
        
        if numMat>1
            partNameSet=[partName,'_',sprintf('%u',matId)];
        else
            partNameSet=partName;
        end
        
        T_elem_sub=cell(n_steps+2,1);
        T_elem_sub(1,1)={['		<Elements mat="',sprintf('%u',matId),'" name="',partNameSet,'" type="',E_type,'">']};
        %               <Elements mat="1" name="Cube_1" type="hex8">
        for q_e2=1:1:n_steps
            %              <elem id="11">   13,   19,   20,   14,   49,   55,   56,   50</elem>
            
            e_ind=elementIndices_mat(q_e2);
            
            T_elem_sub(q_e2+1,1)={['			<elem id="',sprintf('%u',e_ind),'">',sprintf(t_form,E_mat(q_e2,:)),'</elem>']};
            
            if FEB_struct.disp_opt==1 && rem(q_e2,round(n_steps/10))==0
                waitbar(q_e2/n_steps);
            end
        end
        T_elem_sub(end,1)={'		</Elements>'};
        
        T_elem(end+1:end+numel(T_elem_sub),1)=T_elem_sub;
    end
    
    if FEB_struct.disp_opt==1
        close(hw); drawnow;
    end
end



%% MeshData field

% <MeshData>
%       <ElementData var="shell thickness" elem_set="part1">
%           <elem lid="1">1.0, 1.0, 1.0, 1.0</elem>
%           <elem lid="2">1.0, 1.0, 1.0, 1.0</elem>
%       </ElementData >
%       <ElementData var="pre_stretch" elem_set="part2">
%           <elem lid="1">1.05</elem>
%           <elem lid="2">1.05</elem>
%       </ElementData >
% </MeshData>


% ElementData for thickness / mataxis information

T_meshData={};
if isfield(FEB_struct.Geometry,'MeshData')
    disp('----> Adding MeshData field');
    
    if isfield(FEB_struct.Geometry.MeshData,'ElementData')
        disp('----> Adding element data');
        T_elementData={};
        for q_set=1:1:numel(FEB_struct.Geometry.MeshData.ElementData)
            
            T_elementDataMiddle={};
            
            for q_el=1:1:size(FEB_struct.Geometry.MeshData.ElementData{q_set}.Val,1)
                t_form=repmat(' %6.7e,',1,size(FEB_struct.Geometry.MeshData.ElementData{q_set}.Val,2)); t_form=t_form(1:end-1); %text format for sprintf for current element type
                v_text=sprintf(t_form,FEB_struct.Geometry.MeshData.ElementData{q_set}.Val(q_el,:));
                T_elementDataMiddle(q_el,1)={['						<elem lid="',sprintf('%u',q_el),'">',v_text,'</elem>']};
            end
            
            %         numElementsPerSet=cellfun(@(x) size(x,1),FEB_struct.Geometry.Elements);
            %         numElementsPerSetCumsum=cumsum(numElementsPerSet);
            %         elementSize2PerSet=cellfun(@(x) size(x,2),FEB_struct.Geometry.Elements);
            %
            %         elementSize2=zeros(sum(numElementsPerSet),1);
            %         c=1;
            %         for q=1:1:numel(FEB_struct.Geometry.Elements) %element sets
            %             elementSize2(c:numElementsPerSetCumsum(q))=elementSize2PerSet(q);
            %             c=numElementsPerSetCumsum(q)+1;
            %         end
            %
            %
            %         if isfield(FEB_struct.Geometry.ElementData,'Thickness')
            %             elementIndicesThickness=FEB_struct.Geometry.ElementData.IndicesForThickness;
            %             disp('----> Thickness data entries found');
            %         else
            %             elementIndicesThickness=[];
            %         end
            %
            %         if isfield(FEB_struct.Geometry.ElementData,'MatAxis')
            %             elementIndicesMatAxis=FEB_struct.Geometry.ElementData.MatAxis.ElementIndices;
            %             disp('----> MatAxis data entries found');
            %         else
            %             elementIndicesMatAxis=[];
            %         end
            %
            %         if isfield(FEB_struct.Geometry.ElementData,'Fiber')
            %             elementIndicesFiber=FEB_struct.Geometry.ElementData.Fiber.ElementIndices;
            %             disp('----> MatAxis data entries found');
            %         else
            %             elementIndicesFiber=[];
            %         end
            %
            %         elementIndices=unique([elementIndicesThickness(:);elementIndicesMatAxis(:);elementIndicesFiber(:)]);
            %         elementIndices=elementIndices(:)';
            %
            %         T_elementDataMiddle=cell((numel(elementIndices)*2+numel(elementIndicesThickness)+4*numel(elementIndicesMatAxis)+numel(elementIndicesFiber)),1);
            %
            %         n_steps=numel(elementIndices);
            %         if FEB_struct.disp_opt==1
            %             hw = waitbar(0,['Adding ElementData entries....']);
            %         end
            %         c=1;
            %         cw=1;
            %         for q_e=elementIndices
            %
            %             T_elementDataMiddle(c,1)={['         <element id="',sprintf('%u',q_e),'">']};
            %
            %             %Thickness entries
            %             if ismember(q_e,elementIndicesThickness)
            %
            %                 t_form=repmat(' %6.7e,',1,elementSize2(q_e)); t_form=t_form(1:end-1); %text format for sprintf for current element type
            %
            %                 thicknessData=FEB_struct.Geometry.ElementData.Thickness(FEB_struct.Geometry.ElementData.IndicesForThickness==q_e)*ones(1,elementSize2(q_e));
            %                 thicknessEntryText=sprintf(t_form,thicknessData);
            %
            %                 c=c+1;
            %                 T_elementDataMiddle(c,1)={['            <thickness>',thicknessEntryText,'</thickness>']};
            %
            %             end %IF Thickness
            %
            %             %Mat_axis and fibre entries
            %             if ismember(q_e,elementIndicesMatAxis)
            %                 c=c+1;
            %                 T_elementDataMiddle(c,1)={['            <mat_axis>']};
            %
            %                 adData=FEB_struct.Geometry.ElementData.MatAxis.Basis(FEB_struct.Geometry.ElementData.MatAxis.ElementIndices==q_e,:,:);
            %
            %                 c=c+1;
            %                 aData=adData(:,:,1);
            %                 aText=sprintf('%6.7e, %6.7e, %6.7e',aData);
            %                 T_elementDataMiddle(c,1)={['               <a>',aText,'</a>']};
            %
            %                 c=c+1;
            %                 dData=adData(:,:,2);
            %                 dText=sprintf('%6.7e, %6.7e, %6.7e',dData);
            %                 T_elementDataMiddle(c,1)={['               <d>',dText,'</d>']};
            %
            %                 c=c+1;
            %                 T_elementDataMiddle(c,1)={['            </mat_axis>']};
            %
            %             end %IF MatAxis
            %
            %             c=c+1;
            %             T_elementDataMiddle(c,1)={['         </element>']};
            %
            %             if FEB_struct.disp_opt==1 && rem(cw,round(n_steps/10))==0
            %                 waitbar(cw/n_steps);
            %             end
            %             c=c+1;
            %             cw=cw+1;
            %         end %FOR q_e1
            %
            %         if FEB_struct.disp_opt==1
            %             close(hw); drawnow;
            %         end
            %
            
            %Compose elementData cell
            T_elementData(end+1,1)={['				<ElementData var="',FEB_struct.Geometry.MeshData.ElementData{q_set}.Var,'" elem_set="',FEB_struct.Geometry.MeshData.ElementData{q_set}.SetName,'">']};
            T_elementData(end+1:end+numel(T_elementDataMiddle),1)=T_elementDataMiddle;
            T_elementData(end+1,1)={'				</ElementData>'};
        end
    else
        T_elementData={};
    end %IF ElementData
    
    
    if isfield(FEB_struct.Geometry.MeshData,'NodeData')
        disp('----> Adding node data');
        T_nodeData={};
        for q_set=1:1:numel(FEB_struct.Geometry.MeshData.NodeData)
            
            T_nodeDataMiddle={};
            
            for q_el=1:1:size(FEB_struct.Geometry.MeshData.NodeData{q_set}.Val,1)
                t_form=repmat(' %6.7e,',1,size(FEB_struct.Geometry.MeshData.NodeData{q_set}.Val,2)); t_form=t_form(1:end-1); %text format for sprintf for current element type
                v_text=sprintf(t_form,FEB_struct.Geometry.MeshData.NodeData{q_set}.Val(q_el,:));
                T_nodeDataMiddle(q_el,1)={['						<node lid="',sprintf('%u',q_el),'">',v_text,'</node>']};
            end
            
            %Compose nodeData cell
            %             T_nodeData(end+1,1)={['				<NodeData var="',FEB_struct.Geometry.MeshData.NodeData{q_set}.Var,'" node_set="',FEB_struct.Geometry.MeshData.NodeData{q_set}.SetName,'">']};
            T_nodeData(end+1,1)={['				<NodeData name="',FEB_struct.Geometry.MeshData.NodeData{q_set}.Var,'" node_set="',FEB_struct.Geometry.MeshData.NodeData{q_set}.SetName,'">']};
            T_nodeData(end+1:end+numel(T_nodeDataMiddle),1)=T_nodeDataMiddle;
            T_nodeData(end+1,1)={'				</NodeData>'};
        end
    else
        T_nodeData={};
    end %IF NodeData
    
    %Compose meshData cell
    T_meshData(1,1)={'		<MeshData>'};
    if ~isempty(T_elementData)
        T_meshData(end+1:end+numel(T_elementData),1)=T_elementData;
    end
    if ~isempty(T_nodeData)
        T_meshData(end+1:end+numel(T_nodeData),1)=T_nodeData;
    end
    T_meshData(end+1,1)={'		</MeshData>'};
    
else
    T_meshData={};
end

%% Adding surface field

T_surf={};
if isfield(FEB_struct.Geometry,'Surface')
    disp('----> Adding surface field');
    
    for q_set=1:1:numel(FEB_struct.Geometry.Surface)
        
        F=FEB_struct.Geometry.Surface{q_set}.Set; %Faces
        numFaces=size(F,1);
        surfaceType=FEB_struct.Geometry.Surface{q_set}.Type; %Surface type
        
        %Create surface level
        if ~isfield(FEB_struct.Geometry.Surface{q_set},'Name')
            FEB_struct.Geometry.Surface{q_set}.Name=['Surface_',sprintf('%u',q_set)]; %Surface type
        end
        surfaceName=FEB_struct.Geometry.Surface{q_set}.Name; %Surface type
        
        T_surf_sub=cell(numFaces+2,1);
        T_surf_sub(1,1)={['		<Surface name="',surfaceName,'">']}; %      <Surface name="Contact_master">
        
        for q_face=1:1:numFaces %For all faces
            t_form=repmat('   %u,',1,size(F,2)); t_form=t_form(1:end-1); %text format for sprintf for current element type
            T_surf_sub(q_face+1,1)={['         <',surfaceType,' id="',sprintf('%u',q_face),'">',sprintf(t_form,F(q_face,:)),'</',surfaceType,'>']}; %          <tri3 id="1">   4906,   5101,   5100</tri3>
        end
        T_surf_sub(end,1)={'		</Surface>'};
        
        T_surf(end+1:end+numel(T_surf_sub),1)=T_surf_sub;
        
    end
else
    T_surf={};
end

%% Adding NodeSet section

T_nodeSet={};

if isfield(FEB_struct.Geometry,'NodeSet')
    disp('----> Adding NodeSet field');
    
    for q_set=1:1:numel(FEB_struct.Geometry.NodeSet)
        
        nodeSetIndices=FEB_struct.Geometry.NodeSet{q_set}.Set; %Node indices
        nodeSetIndices=nodeSetIndices(:)'; %Row vector
        numNodes=numel(nodeSetIndices);
        
        if ~isfield(FEB_struct.Geometry.NodeSet{q_set},'Name')
            FEB_struct.Geometry.NodeSet{q_set}.Name=['NodeSet_',sprintf('%u',q_set)]; %Surface type
        end
        nodeSetName=FEB_struct.Geometry.NodeSet{q_set}.Name; %Node set name
        
        T_nodeSet_sub=cell(numNodes+2,1);
        T_nodeSet_sub(1,1)={['		<NodeSet name="',nodeSetName,'">']}; %      <Surface name="Contact_master">
        
        for q_node=1:1:numNodes %For all faces
            T_nodeSet_sub(q_node+1,1)={['         <node id="',sprintf('%u',nodeSetIndices(q_node)),'"/>']}; % <node id="1"/>
        end
        T_nodeSet_sub(end,1)={'		</NodeSet>'};
        T_nodeSet(end+1:end+numel(T_nodeSet_sub),1)=T_nodeSet_sub;
        
    end
    
end

%% Adding ElementSet section

T_elementSet={};

if isfield(FEB_struct.Geometry,'ElementSet')
    disp('----> Adding ElementSet field');
    
    for q_set=1:1:numel(FEB_struct.Geometry.ElementSet)
        
        elementSetIndices=FEB_struct.Geometry.ElementSet{q_set}.Set; %Node indices
        elementSetIndices=elementSetIndices(:)'; %Row vector
        numNodes=numel(elementSetIndices);
        
        if ~isfield(FEB_struct.Geometry.ElementSet{q_set},'Name')
            FEB_struct.Geometry.ElementSet{q_set}.Name=['ElementSet_',sprintf('%u',q_set)]; %Surface type
        end
        elementSetName=FEB_struct.Geometry.ElementSet{q_set}.Name; %Node set name
        
        T_elementSet_sub=cell(numNodes+2,1);
        T_elementSet_sub(1,1)={['		<ElementSet name="',elementSetName,'">']}; %      <Surface name="Contact_master">
        
        for q_node=1:1:numNodes %For all faces
            T_elementSet_sub(q_node+1,1)={['         <elem id="',sprintf('%u',elementSetIndices(q_node)),'"/>']}; % <node id="1"/>
        end
        T_elementSet_sub(end,1)={'		</ElementSet>'};
        T_elementSet(end+1:end+numel(T_elementSet_sub),1)=T_elementSet_sub;
        
    end
    
end

%% Adding SurfacePair section

T_surfacePair={};

if isfield(FEB_struct.Geometry,'SurfacePair')
    disp('----> Adding SurfacePair field');
    
    for q_set=1:1:numel(FEB_struct.Geometry.SurfacePair)
        
        if ~isfield(FEB_struct.Geometry.SurfacePair{q_set},'Name')
            FEB_struct.Geometry.SurfacePair{q_set}.Name=['SurfacePair_',sprintf('%u',q_set)]; %Surface type
        end
        surfacePairName=FEB_struct.Geometry.SurfacePair{q_set}.Name; %Node set name
        
        T_surfacePair_sub=cell(4,1);
        T_surfacePair_sub(1,1)={['		<SurfacePair name="',surfacePairName,'">']}; %      <Surface name="Contact_master">
        T_surfacePair_sub(2,1)={['         <master surface="',FEB_struct.Geometry.SurfacePair{q_set}.Master,'"/>']}; % <node id="1"/>
        T_surfacePair_sub(3,1)={['         <slave surface="',FEB_struct.Geometry.SurfacePair{q_set}.Slave,'"/>']}; % <node id="1"/>
        T_surfacePair_sub(end,1)={'		</SurfacePair>'};
        T_surfacePair(end+1:end+numel(T_surfacePair_sub),1)=T_surfacePair_sub;
        
    end
    
end

%% Compose text cell
totalTextCell=T_top;
totalTextCell(end+1:end+numel(T_node))=T_node;
totalTextCell(end+1:end+numel(T_elem))=T_elem;

if ~isempty(T_meshData)
    totalTextCell(end+1:end+numel(T_meshData))=T_meshData;
end

if ~isempty(T_surf)
    totalTextCell(end+1:end+numel(T_surf))=T_surf;
end

if ~isempty(T_nodeSet)
    totalTextCell(end+1:end+numel(T_nodeSet))=T_nodeSet;
end

if ~isempty(T_elementSet)
    totalTextCell(end+1:end+numel(T_elementSet))=T_elementSet;
end

if ~isempty(T_surfacePair)
    totalTextCell(end+1:end+numel(T_surfacePair))=T_surfacePair;
end

totalTextCell(end+1:end+numel(T_end))=T_end;

%% Export text cell to .feb file
cell2txtfile(FEB_struct.run_filename,totalTextCell,1); %Export to text file

%% Reimport XML type
domNode = xmlread(FEB_struct.run_filename); %Reimport docNode

end
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2017  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
