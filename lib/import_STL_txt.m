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
