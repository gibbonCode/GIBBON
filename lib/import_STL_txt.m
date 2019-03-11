function [stlStruct] = import_STL_txt(fileName)

% function [stlStruct] = import_STL_txt(fileName)
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
% 2016/02/24 Updated commenting and header
%------------------------------------------------------------------------

T=txtfile2cell(fileName);

logicSolids= ~cellfun(@isempty,strfind(T,'solid'));
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
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
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
