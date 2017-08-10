function docNode=add_fiber_dir_FEB(docNode,FEB_struct,disp_opt)

%This function adds fiber direction data to a FEBIo .feb input file
%
% The following forms are supported:
%
% For opt='transiso':
%<ElementData>
%   <element id="E_ind(i)">
%         <fiber>Vf(i,1),Vf(i,2),Vf(i,3)</fiber>
%   </element>
%</ElementData>
%
%
% For opt='ortho':
%<ElementData>
%   <element id="E_ind(i)">
%      <mat_axis>
%             <a>Vf(i,1,1),Vf(i,2,1),Vf(i,3,1)</a>
%             <d>Vf(i,1,2),Vf(i,2,2),Vf(i,3,2)</d>
%       </mat_axis>
%   </element>
%</ElementData>
%
% 25/01/2012, Kevin Mattheus Moerman
% kevinmoerman@gmail.com

%%

warning('add_fiber_dir_FEB will be replaced by addMatAxisFibreElementData_FEB in future releases');

docNode=addMatAxisFibreElementData_FEB(docNode,FEB_struct);


% Vf=FEB_struct.Geometry.ElementData.MatAxis.Basis;
% E_ind=FEB_struct.Geometry.ElementData.MatAxis.ElementIndices;
% 
% %Get Geometry level
% GEONode = docNode.getElementsByTagName('Geometry').item(0);
% docRootNode =  GEONode;
% 
% %Check for ElementData level
% if GEONode.getElementsByTagName('ElementData').getLength==0; %ElementData level does not exist yet
%     %Adding ElementData level
%     parent_node = docNode.createElement('ElementData');
%     parent_node = docRootNode.appendChild(parent_node);
%     new_entry=1;
% else %ElementData level already exists
%     %Finding existing element entries
%     no_entries=docRootNode.getElementsByTagName('element').getLength;
%     ID_no=zeros(1,no_entries);
%     for i=0:1:no_entries-1
%         ID_no(i+1)=str2double(docRootNode.getElementsByTagName('element').item(i).getAttribute('id').toCharArray()');
%     end
%     new_entry=0;
% end
% docRootNode =  GEONode.getElementsByTagName('ElementData').item(0);
% 
% %Creating ElementData entries
% if disp_opt==1;
%     hw = waitbar(0,'Creating MatAxis entries...');
% end
% disp('----> Creating MatAxis entries')
% 
% n_steps = numel(E_ind);
% for i=1:n_steps
%     if new_entry==1 || ~any(ID_no==E_ind(i)); %Need to add element level
%         element_node = docNode.createElement('element');
%         element_node = docRootNode.appendChild(element_node);
%         attr = docNode.createAttribute('id');
%         attr.setNodeValue(num2str(E_ind(i)));
%         element_node.setAttributeNode(attr);
%     else %element level already exists need to check for fiber entries
%         element_node=docRootNode.getElementsByTagName('element').item(find(ID_no==E_ind(i))-1);
%     end
%     
%     if size(Vf,3)==1
%         %case 'transiso'
%             v_text=sprintf('%6.7e, %6.7e, %6.7e',Vf(i,:));
%             if element_node.getElementsByTagName('fiber').getLength==0
%                 %Adding fiber level
%                 fiber_node = docNode.createElement('fiber');
%                 element_node.appendChild(fiber_node);
%                 %Adding fiber data text
%                 fiber_node.appendChild(docNode.createTextNode(v_text));
%             else %Fiber entries exist, overwriting existing
%                 element_node.getElementsByTagName('fiber').item(0).getFirstChild.setData(char(v_text));
%             end
%     else
%         %case 'ortho'
%             a_text=sprintf('%6.7e, %6.7e, %6.7e',Vf(i,:,1));
%             d_text=sprintf('%6.7e, %6.7e, %6.7e',Vf(i,:,2));
%             
%             if element_node.getElementsByTagName('mat_axis').getLength==0
%                 %Adding mat_axis level
%                 mat_axis_node = docNode.createElement('mat_axis');
%                 element_node.appendChild(mat_axis_node);
%                 %Adding a level
%                 a_node = docNode.createElement('a');
%                 mat_axis_node.appendChild(a_node);
%                 %Adding a data text
%                 a_node.appendChild(docNode.createTextNode(a_text));
%                 %Adding d level
%                 d_node = docNode.createElement('d');
%                 mat_axis_node.appendChild(d_node);
%                 %Adding d data text
%                 d_node.appendChild(docNode.createTextNode(d_text));
%             else %Mat_axis entries exist, overwriting existing
%                 %Finding mat_axis level
%                 mat_axis_node=element_node.getElementsByTagName('mat_axis').item(0);
%                 %Overwriting a level
%                 mat_axis_node.getElementsByTagName('a').item(0).getFirstChild.setData(char(a_text));
%                 %Overwriting d level
%                 mat_axis_node.getElementsByTagName('d').item(0).getFirstChild.setData(char(d_text));
%             end
%     end
%     if disp_opt==1;
%         waitbar(i/n_steps);
%     end
% end
% 
% if disp_opt==1;
%     close(hw);
% end
% 
% end
 
%% <-- GIBBON footer text --> 
% 
%     GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
%     image segmentation, image-based modeling, meshing, and finite element
%     analysis.
% 
%     Copyright (C) 2017  Kevin Mattheus Moerman
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
