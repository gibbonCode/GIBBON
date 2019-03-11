function [meshStruct]=hexMeshSphere(optionStruct)

%% Parse input

%Compare to default structure
defaultOptionStruct.sphereRadius=1;
defaultOptionStruct.coreRadius=defaultOptionStruct.sphereRadius/2;
defaultOptionStruct.numElementsMantel=3;
defaultOptionStruct.numElementsCore=3; 
defaultOptionStruct.makeHollow=0;
defaultOptionStruct.outputStructType=1;
defaultOptionStruct.cParSmooth.Method='LAP';
defaultOptionStruct.cParSmooth.LambdaSmooth=0.5;
defaultOptionStruct.cParSmooth.n=5;
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1);

% Get input parameters
sphereRadius=optionStruct.sphereRadius;
coreRadius=optionStruct.coreRadius;
numElementsMantel=optionStruct.numElementsMantel;
numElementsCore=optionStruct.numElementsCore;
makeHollow=optionStruct.makeHollow;
outputStructType=optionStruct.outputStructType;

%% CREATING A HEXAHEDRAL MESH CUBE

%Create box 1
sphereDim=1/sqrt(3)*2*coreRadius*ones(1,3); 
sphereEl=numElementsCore*ones(1,3); %Number of elements
[boxMeshStruct]=hexMeshBox(sphereDim,sphereEl);
E_core=boxMeshStruct.E;
V_core=boxMeshStruct.V;
Fb=boxMeshStruct.Fb;
faceBoundaryMarkerBox=boxMeshStruct.faceBoundaryMarker;

indBoundary=unique(Fb(:));
V_core_boundary=V_core(indBoundary,:);

%% MAPPING OUTER SURFACE TO A SPHERE

[azimuth,elevation,r] = cart2sph(V_core_boundary(:,1),V_core_boundary(:,2),V_core_boundary(:,3));
[V_core_boundary(:,1),V_core_boundary(:,2),V_core_boundary(:,3)] = sph2cart(azimuth,elevation,coreRadius.*ones(size(r)));
V_core(indBoundary,:)=V_core_boundary;

%% Adding mantel

mantelThickness=sphereRadius-coreRadius;
[Fq,Vq,~]=patchCleanUnused(Fb,V_core);
[E_mantel,V_mantel,F_mantel_inner,F_mantel_outer]=quadThick(Fq,Vq,1,mantelThickness,numElementsMantel);

%Fix outer radii
indBoundary=unique(F_mantel_outer(:));
V_mantel_boundary=V_mantel(indBoundary,:);

[azimuth,elevation,r] = cart2sph(V_mantel_boundary(:,1),V_mantel_boundary(:,2),V_mantel_boundary(:,3));
[V_mantel_boundary(:,1),V_mantel_boundary(:,2),V_mantel_boundary(:,3)] = sph2cart(azimuth,elevation,sphereRadius.*ones(size(r)));
V_mantel(indBoundary,:)=V_mantel_boundary;

%%
if makeHollow==1
    ET=E_mantel; 
    VT=V_mantel; 
    FTb=[F_mantel_outer; F_mantel_inner];
    faceBoundaryMarker=[ones(size(F_mantel_outer,1),1); 2*ones(size(F_mantel_inner,1),1)];
elseif makeHollow==0    
    % Merging node sets
    VT=[V_core;V_mantel];
    ET=[E_core;E_mantel+size(V_core,1)];
    [FT]=element2patch(ET,[],'hex8');
    [FT,VT,~,ind2]=mergeVertices(FT,VT);    
    ET=ind2(ET);
    F_mantel_outer=ind2(F_mantel_outer+size(V_core,1));
%     F_mantel_inner=ind2(F_mantel_inner+size(V_core,1));
    FTb=F_mantel_outer;
    faceBoundaryMarker=ones(size(F_mantel_outer,1),1);
end

%Get faces
[FT,~]=element2patch(ET,[],'hex8');

%% Smoothing

optionStruct.cParSmooth.RigidConstraints=unique(FTb(:));

[F,~,~]=uniqueIntegerRow(FT);

if optionStruct.cParSmooth.n>0
    [VT]=tesSmooth(F,VT,[],optionStruct.cParSmooth);
end

%% Collect output

switch outputStructType
    case 1
        meshStruct.E=ET;
        if makeHollow==0
            meshStruct.elementRegionLabel=[1*ones(size(E_core,1),1); 2*ones(size(E_mantel,1),1);];
        end
        meshStruct.V=VT;
        meshStruct.F=FT;
        meshStruct.Fb=FTb;
        meshStruct.faceBoundaryMarker=faceBoundaryMarker;
        meshStruct.faceBoundaryMarkerBox=faceBoundaryMarkerBox;
    case 2
        meshStruct.nodes=VT;
        meshStruct.facesBoundary=FTb;
        meshStruct.boundaryMarker=faceBoundaryMarker;
        meshStruct.faces=FT;
        meshStruct.elements=ET;        
        if makeHollow==0
            meshStruct.elementMaterialID=[1*ones(size(E_core,1),1); 2*ones(size(E_mantel,1),1);];
        else
            meshStruct.elementMaterialID=ones(size(ET,1),1);
        end
        meshStruct.faceMaterialID=ones(size(meshStruct.faces,1),1);
end


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
