function [meshStruct]=hexMeshBox(varargin)

% function [meshStruct]=hexMeshBox(boxDim,boxEl,outputStructType)
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
% 2018/01/23 Added option to export same mesh structure as tetgen so it is
% compatible with the meshView function
%------------------------------------------------------------------------

%% Parse input 

switch nargin
    case 1
        boxDim=varargin{1};
        boxEl=[10 10 10];
        outputStructType=1;
    case 2
        boxDim=varargin{1};
        boxEl=varargin{2};
        outputStructType=1;
    case 3
        boxDim=varargin{1}; 
        boxEl=varargin{2}; 
        outputStructType=varargin{3}; 
end

%%
dX=boxDim(1); 
dY=boxDim(2);
dZ=boxDim(3);
nX=boxEl(1); 
nY=boxEl(2);
nZ=boxEl(3);

%Creating 2D slice
[X2,Y2]=meshgrid(linspace(-dX/2,dX/2,nX+1),linspace(-dY/2,dY/2,nY+1)); %mesh of single slice
[F2,V2,~] = grid2patch(X2,Y2,zeros(size(X2)),zeros(size(X2))); %Convert to patch data (quadrilateral faces)

%Offsetting to create 3D model
z_range=linspace(-dZ/2,dZ/2,nZ+1); %Offset steps
V=repmat(V2,numel(z_range),1); %replicating coordinates
Z_add=ones(size(V2,1),1)*z_range;
Z_add=Z_add(:);
V(:,3)=V(:,3)+Z_add; %Altering z coordinates

%Creating element definitions
q=(((1:(nZ))-1).*size(V2,1));
Q=q(ones(size(F2,1),1),:);
Qc=Q(:);
QQc=Qc(:,ones(1,size(F2,2)));

q=((1:(nZ)).*size(V2,1));
Q=q(ones(size(F2,1),1),:);
Qc=Q(:);
QQc2=Qc(:,ones(1,size(F2,2)));

E=[repmat(F2,[nZ 1])+(QQc) repmat(F2,[nZ 1])+(QQc2); ];

%Convert elements to faces
[F,~]=element2patch(E,[],'hex8');

%Find boundary faces
[indFree]=freeBoundaryPatch(F);
Fb=F(indFree,:);

%Create faceBoundaryMarkers based on normals
[N]=patchNormal(Fb,V); %N.B. Change of convention changes meaning of front, top etc.

faceBoundaryMarker=zeros(size(Fb,1),1);

faceBoundaryMarker(N(:,1)<-0.5)=1; %Left
faceBoundaryMarker(N(:,1)>0.5)=2; %Right
faceBoundaryMarker(N(:,2)<-0.5)=3; %Front
faceBoundaryMarker(N(:,2)>0.5)=4; %Back
faceBoundaryMarker(N(:,3)<-0.5)=5; %Bottom
faceBoundaryMarker(N(:,3)>0.5)=6; %Top

%% Collect output

switch outputStructType
    case 1
        meshStruct.E=E;
        meshStruct.V=V;
        meshStruct.F=F;
        meshStruct.indFree=indFree;
        meshStruct.Fb=Fb;
        meshStruct.faceBoundaryMarker=faceBoundaryMarker;
    case 2
        meshStruct.nodes=V;
        meshStruct.facesBoundary=Fb;
        meshStruct.boundaryMarker=faceBoundaryMarker;
        meshStruct.faces=F;
        meshStruct.elements=E;
        meshStruct.elementMaterialID=ones(size(E,1),1);
        meshStruct.faceMaterialID=ones(size(meshStruct.faces,1),1);
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
