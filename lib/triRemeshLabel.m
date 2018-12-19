function [Fb,Vb,Cb]=triRemeshLabel(F,V,pointSpacing)

% function [Fb,Vb,Cb]=triRemeshLabel(F,V,pointSpacing)
% ------------------------------------------------------------------------
% This function uses the PDE toolbox remeshing tools to remesh and label
% the input surfaces. 
%
%
% ------------------------------------------------------------------------

%%

% Initialize model
modelTemp = createpde('thermal','steadystate'); 

% Assign geometry
modelTemp.Geometry = geometryFromMesh(modelTemp,V',F');

% Mesh
modelMesh=generateMesh(modelTemp,'Hmax',pointSpacing,'GeometricOrder','linear');

% Access PDE-toolbox created mesh
V=modelTemp.Mesh.Nodes';
E=modelTemp.Mesh.Elements';

% Get mesh faces
F=element2patch(E,[],'tet4');

% Get boundary faces
indBoundary=tesBoundary(F,V);
Fb=F(indBoundary,:);
Vb=V;

% Get face labels
Cb=zeros(size(Fb,1),1);
for q=1:1:modelTemp.Geometry.NumFaces
     indFind = findNodes(modelMesh,'region','Face',q);     
     Cb(all(ismember(Fb,indFind),2))=q;     
end

% Remove unused nodes
[Fb,Vb]=patchCleanUnused(Fb,Vb);

% Remove three connected triangles
[Fb,Vb,Cb]=triSurfRemoveThreeConnect(Fb,Vb,Cb);

