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
modelMesh=generateMesh(modelTemp,'Hmin',pointSpacing*0.9,'Hmax',pointSpacing*1.1,'GeometricOrder','linear');

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
% [Fb,Vb,Cb]=triSurfRemoveThreeConnect(Fb,Vb,Cb);

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
