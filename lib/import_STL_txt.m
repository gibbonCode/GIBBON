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
