function [stlStruct] = import_STL(fileName)

% function [stlStruct] = import_STL_txt(fileName)
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2016/02/24
% To do: create proper test to check for binary files
%------------------------------------------------------------------------

%% Read text data

T=txtfile2cell(fileName);
logicSolids= ~cellfun(@isempty,strfind(T,'solid'));

try
    if any(logicSolids)
        [stlStruct]=parse_STL_txt(T,logicSolids);
    else
        [stlStruct]=parse_STL_bin(fileName);
    end
catch
    [stlStruct]=parse_STL_bin(fileName);
end

end

%%

function [stlStruct]=parse_STL_txt(T,logicSolids)

numSolids=nnz(logicSolids)/2;
stlStruct.solidNames=cell(1,numSolids);
stlStruct.solidVertices=cell(1,numSolids);
stlStruct.solidFaces=cell(1,numSolids);
stlStruct.solidNormals=cell(1,numSolids);

for q=1:1:numSolids
    
    %Find current solid
    indStart=find(logicSolids,1);
    logicSolids(indStart)=0;
    indEnd=find(logicSolids,1);
    logicSolids(indEnd)=0;
    
    T_solid=T(indStart:indEnd); %Current solid text
    
    solidName = sscanf(T_solid{1}, 'solid %s'); %Solid name
    
    %Get vertices
    logicVertices= ~cellfun(@isempty,strfind(T_solid,'vertex'));
    numVertices=nnz(logicVertices);
    T_vertices=T_solid(logicVertices);
    V=cell2mat(cellfun(@(x) sscanf(x,'    vertex %f %f %f')',T_vertices,'UniformOutput',0));
    
    %Get face normals
    logicNormals= ~cellfun(@isempty,strfind(T_solid,'facet normal'));
    T_normals=T_solid(logicNormals);
    N=cell2mat(cellfun(@(x) sscanf(x,'  facet normal %f %f %f')',T_normals,'UniformOutput',0));
    
    %Create F
    F=reshape(1:numVertices,[3 numVertices/3])';
    
    stlStruct.solidNames{q}=solidName;
    stlStruct.solidVertices{q}=V;
    stlStruct.solidFaces{q}=F;
    stlStruct.solidNormals{q}=N;
    
end
end

function [stlStruct]=parse_STL_bin(fileName)

% See also  import_STL_bin

%% Access file

fid=fopen(fileName, 'r'); %Open the file, assumes STL Binary format.
if fid == -1
    error('File could not be opened, check name or path.')
end

%% Read title
solidNameBin=fread(fid,80,'uchar=>schar'); % Read file title
solidName = char(solidNameBin');

%% Get number of faces
numFaces=fread(fid,1,'int32'); % Read number of Faces

%% Read remainder

T = fread(fid,inf,'uint8=>uint8'); % read the remaining values
fclose(fid); %Close file

%% Process import

% Each facet is 50 bytes
%  - Three single precision values specifying the face normal vector
%  - Three single precision values specifying the first vertex (XYZ)
%  - Three single precision values specifying the second vertex (XYZ)
%  - Three single precision values specifying the third vertex (XYZ)
%  - Two color bytes (possibly zeroed)

% 3 dimensions x 4 bytes x 4 vertices = 48 bytes for triangle vertices
% 2 bytes = color (if color is specified)

trilist = 1:48;

ind = reshape(repmat(50*(0:(numFaces-1)),[48,1]),[1,48*numFaces])+repmat(trilist,[1,numFaces]);
coordNormData = reshape(typecast(T(ind),'single'),[3,4,numFaces]);

N=squeeze(coordNormData(:,1,:))';
N=double(N);

V=coordNormData(:,2:4,:);
V=reshape(V,[3,3*numFaces]);
V=double(V)';

F=reshape(1:3*numFaces,[3,numFaces])';


c0 = typecast(T(49:50),'uint16');
if (bitget(c0(1),16)==1)
    trilist = 49:50;
    ind = reshape(repmat(50*(0:(numFaces-1)),[2,1]),[1,2*numFaces])+repmat(trilist,[1,numFaces]);
    c0 = reshape(typecast(T(ind),'uint16'),[1,numFaces]);
    
    r=bitshift(bitand(2^16-1, c0),-10);
    g=bitshift(bitand(2^11-1, c0),-5);
    b=bitand(2^6-1, c0);
    C=[r; g; b]';
else
    C = zeros(numFaces,3);
end


%% Arrange output in structure array
stlStruct.solidNames{1}=solidName;
stlStruct.solidVertices{1}=V;
stlStruct.solidFaces{1}=F;
stlStruct.solidNormals{1}=N;
stlStruct.solidColors{1}=C;

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
