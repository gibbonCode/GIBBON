function [meshStruct]=hexMeshCylinder(varargin)

% function [meshStruct]=hexMeshCylinder(cylRadius,cylLength,pointSpacing)
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2021/02/26 Created
% 2021/05/08 @Fireedman (Luis Antonio Aguilar) added "end" to function
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 0
        cylRadius=[];
        cylLength=[];
        pointSpacing=[];
    case 1
        cylRadius=varargin{1};
        cylLength=[];
        pointSpacing=[];
    case 2
        cylRadius=varargin{1};
        cylLength=varargin{2};
        pointSpacing=[];
    case 3
        cylRadius=varargin{1};
        cylLength=varargin{2};
        pointSpacing=varargin{3};
end

%Check for empty variables
if isempty(cylRadius)
    cylRadius=1; %Unit radius
end

if isempty(cylLength)
    cylLength=cylRadius; %Height the same as radius
end

if isempty(pointSpacing)
    pointSpacing=(cylRadius*2*pi)/10; %One tenth of circumference
end

if numel(pointSpacing)==1
    pointSpacing=pointSpacing.*ones(1,2);
end

%%
%Set number of elements for core
numElementsCore=ceil(((cylRadius*2*pi)./pointSpacing(1))/4);

if numElementsCore<1
    numElementsCore=1;
end

%Set number of nodes allong circumference
numElementsHeight=ceil(cylLength./pointSpacing(2));
if numElementsHeight<1
    numElementsHeight=1;
end

%% Create quad mesh for disc

%Raw quad mesh
[Fs,Vs]=discQuadMesh(numElementsCore,cylRadius,0.6);

%Smoothen
Eb=patchBoundary(Fs);
controlParSmooth.n=25;
controlParSmooth.Method='LAP';
controlParSmooth.RigidConstraints=unique(Eb(:));
[Vs]=patchSmooth(Fs,Vs,[],controlParSmooth);

%% Thicken to elements

%Thicken
[E,V,Fp1,Fp2]=patchThick(Fs,Vs,1,cylLength,numElementsHeight);
V(:,3)=V(:,3)-cylLength/2;

%Get element faces
[F,~]=element2patch(E,[],'hex8');

%Find boundary faces
[indFree]=freeBoundaryPatch(F);
Fb=F(indFree,:);

%Assign boundary labels (or colors)
faceBoundaryMarker=zeros(size(Fb,1),1); %Side of cylinder
faceBoundaryMarker(all(ismember(Fb,Fp1),2))=1; %Original (e.g. top)
faceBoundaryMarker(all(ismember(Fb,Fp2),2))=2; %Thickened (e.g. bottom)

%% Collect output

meshStruct.nodes=V;
meshStruct.facesBoundary=Fb;
meshStruct.boundaryMarker=faceBoundaryMarker;
meshStruct.faces=F;
meshStruct.elements=E;
meshStruct.elementMaterialID=ones(size(E,1),1);
meshStruct.faceMaterialID=ones(size(F,1),1);

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
% Copyright (C) 2006-2021 Kevin Mattheus Moerman and the GIBBON contributors
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
