function domNode=addMatAxisFibreElementData_FEB(domNode,FEB_struct)

% function domNode=addMatAxisFibreElementData_FEB(domNode,FEB_struct)
% ------------------------------------------------------------------------
%
%
% 2016/05/12 Taken out of addGeometryLevel_FEB as seperate function to
% replace add_fiber_dir_FEB
% 2016/09/09 Added error if no geometry section exists
% TO DO fix support for fiber and other types (currently only mat_axis
% supported)
%------------------------------------------------------------------------

%%

if ~isfield(FEB_struct,'disp_opt')
    FEB_struct.disp_opt=0;
end
    
%%

Vf=FEB_struct.Geometry.ElementData.MatAxis.Basis;
E_ind=FEB_struct.Geometry.ElementData.MatAxis.ElementIndices;

%Get Geometry level
GEONode = domNode.getElementsByTagName('Geometry').item(0);

if ~isempty(GEONode)
    
    docRootNode =  GEONode;
    
    %Check for ElementData level
    if GEONode.getElementsByTagName('ElementData').getLength==0 %ElementData level does not exist yet
        %Adding ElementData level
        parent_node = domNode.createElement('ElementData');
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
    if FEB_struct.disp_opt==1
        hw = waitbar(0,'Creating MatAxis entries...');
    end
    disp('----> Creating MatAxis entries')
    
    n_steps = numel(E_ind);
    for i=1:n_steps
        if new_entry==1 || ~any(ID_no==E_ind(i)) %Need to add element level
            element_node = domNode.createElement('element');
            element_node = docRootNode.appendChild(element_node);
            attr = domNode.createAttribute('id');
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
                fiber_node = domNode.createElement('fiber');
                element_node.appendChild(fiber_node);
                %Adding fiber data text
                fiber_node.appendChild(domNode.createTextNode(v_text));
            else %Fiber entries exist, overwriting existing
                element_node.getElementsByTagName('fiber').item(0).getFirstChild.setData(char(v_text));
            end
        else
            %case 'ortho'
            a_text=sprintf('%6.7e, %6.7e, %6.7e',Vf(i,:,1));
            d_text=sprintf('%6.7e, %6.7e, %6.7e',Vf(i,:,2));
            
            if element_node.getElementsByTagName('mat_axis').getLength==0
                %Adding mat_axis level
                mat_axis_node = domNode.createElement('mat_axis');
                element_node.appendChild(mat_axis_node);
                %Adding a level
                a_node = domNode.createElement('a');
                mat_axis_node.appendChild(a_node);
                %Adding a data text
                a_node.appendChild(domNode.createTextNode(a_text));
                %Adding d level
                d_node = domNode.createElement('d');
                mat_axis_node.appendChild(d_node);
                %Adding d data text
                d_node.appendChild(domNode.createTextNode(d_text));
            else %Mat_axis entries exist, overwriting existing
                %Finding mat_axis level
                mat_axis_node=element_node.getElementsByTagName('mat_axis').item(0);
                %Overwriting a level
                mat_axis_node.getElementsByTagName('a').item(0).getFirstChild.setData(char(a_text));
                %Overwriting d level
                mat_axis_node.getElementsByTagName('d').item(0).getFirstChild.setData(char(d_text));
            end
        end
        if FEB_struct.disp_opt==1 && rem(i,round(n_steps/10))==0
            waitbar(i/n_steps);
        end
    end
    
    if FEB_struct.disp_opt==1
        close(hw); drawnow;
    end
else
    error('The geometry node appears to be missing');
end
 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
