function export_mesh(fileName,meshStruct)

% function export_mesh(fileName,meshStruct)
% ------------------------------------------------------------------------
% This function exports a .MESH file based on the meshStruct which should
% contain the vertices and elements. 
%
% Elements can be specified as either: 
% meshStruct.elements or meshStruct.E
%
% Vertices can be specified as either: 
% meshStruct.vertices or meshStruct.V or meshStruct.nodes
%
% This function currently only supporte tetrahedral and hexahedral element
% exporting. 
%
% The element type is either inferred from the size of the element matrix
% or can be set using: 
% meshStruct.elementType
% Use meshStruct.elementType='Hexahedral', for 8-noded hexahedral elements
% Use meshStruct.elementType='Tetrahedral', for 4-noded tetrahedral elements
%
% This implementation was created using examples of .MESH files contained
% in: https://github.com/cnr-isti-vclab/HexaLab
% 
%
% Change log: 
% 2021/07/20 Kevin M. Moerman Created for use with https://www.hexalab.net/
%
% To do: 
% - Proper treatment of vertex/element "attributes"
% - Support more .MESH options (surface meshes, FEA types,....)
% - Test for tetrahedral meshes
% - Expand and test for surface meshes
% ------------------------------------------------------------------------

%% Parse input

%Get elements
if isfield(meshStruct,'E')
    E=meshStruct.E;
elseif isfield(meshStruct,'elements')
    E=meshStruct.elements;
end

%Get vertices
if isfield(meshStruct,'V')
    V=meshStruct.V;
elseif isfield(meshStruct,'nodes')
    V=meshStruct.nodes;
elseif isfield(meshStruct,'vertices')
    V=meshStruct.vertices;
end

%Get/set element type
if isfield(meshStruct,'elementType')
    elementType=meshStruct.elementType;
else
    switch size(E,2)
        case 4 % 4 noded tetrahedral elements
            elementType='Tetrahedra';
        case 8 % 8 noded tetrahedral elements
            elementType='Hexahedra';
    end
end

%Replace "FEBio-style" element type designations by .MESH types
if strcmp(elementType,'tet4')
    elementType='Tetrahedra';
elseif strcmp(elementType,'hex8')
    elementType='Hexahedra';
end

%% Open file for writing

fileID = fopen(fileName,'w'); %Open file for writing

%% Write .MESH file version line
% Example: 
%
% MeshVersionFormatted 1

fprintf(fileID,'%s \n','MeshVersionFormatted 1');

%% Write dimensionality line
% Example: 
%
% Dimension 3

fprintf(fileID,'Dimension %d \n',size(V,2)); 

%% Write line specifying number of vertices 
% Example: 
%
% Vertices 12

fprintf(fileID,'Vertices %d \n',size(V,1));

%% Write vertex data field
% Example: 
%
% 0 0 0 1
% -1.0000000000000000e+00 -1.0000000000000000e+00 -1.0000000000000000e+00 1

%Prepare data
V_id=ones(size(V,1),1);
V_field=[V V_id(:)];

%Write field
tForm=repmat('%0.16e ',1,size(V,2)); %Replicated text format
tForm=[tForm(1:end-1),' %i \n']; %Text format with attribute and end of line
fprintf(fileID,tForm,V_field');

%% Write line specifying number of elements and element type
% Example: 
%
% Hexahedra 50

fprintf(fileID,[elementType,' %d \n'],size(E,1));

%% Write element data field

%Prepare data
E_id=ones(size(E,1),1); %1:1:size(E,1);
E_field=[E E_id(:)];

%Write field
tForm=repmat('%d ',1,size(E,2)); %Replicated text format
tForm=[tForm(1:end-1),' %i \n']; %Text format with attribute and end of line
fprintf(fileID,tForm,E_field');

%% Write "End" line
fprintf(fileID,'%s \n','End');

%% Close file for writing
fclose(fileID);

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
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
