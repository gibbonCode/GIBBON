function [varargout]=discQuadMesh(varargin)

% function [F,V,C,indEdge]=discQuadMesh(nElements,r,f)
% ------------------------------------------------------------------------
% This function meshes a circle using quadrilaterial elements. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 10/05/2016 Updated for GIBBON
%------------------------------------------------------------------------

%% Parse input 
switch nargin
    case 1
        nElements=varargin{1};
        r=1;
        f=0.5;
    case 2
        nElements=varargin{1};
        r=varargin{2};
        f=0.5;
    case 3
        nElements=varargin{1};
        r=varargin{2};
        f=varargin{3};
end

%% Creating central regular quad mesh

nElements=nElements+~iseven(nElements);%Force even

[X_centralMesh,Y_centralMesh]=meshgrid(linspace(-1,1,nElements+1));
[F_centralMesh,V_centralMesh] = surf2patch(X_centralMesh,Y_centralMesh,zeros(size(X_centralMesh)));
V_centralMesh=V_centralMesh(:,1:2);

%Edge of central mesh
logicCentralMeshEdge=(X_centralMesh==1)|(Y_centralMesh==1)|(X_centralMesh==-1)|(Y_centralMesh==-1);
nEdge=(nElements*4);

% Scaling radius
[ThetaMesh,RadiusMesh]=cart2pol(V_centralMesh(:,1),V_centralMesh(:,2));
RadiusMesh=f*(1/2)*sqrt(2)*RadiusMesh;
[V_centralMesh(:,1),V_centralMesh(:,2)]=pol2cart(ThetaMesh,RadiusMesh);

%% Creating outer mesh

RadiusOuterEdge=ones(1,nEdge);
ThetaOuterEdge=linspace(0,pi*2,nEdge+1); 
ThetaOuterEdge=ThetaOuterEdge(2:end)-pi;

[xOuterEdge,yOuterEdge]=pol2cart(ThetaOuterEdge,RadiusOuterEdge);
V_outerEdge=[xOuterEdge(:) yOuterEdge(:)];

V_innerEdge=V_centralMesh(logicCentralMeshEdge,:);
[ThetaEdge,RadiusEdge]=cart2pol(V_innerEdge(:,1),V_innerEdge(:,2));
[ThetaEdge,sortInd]=sort(ThetaEdge);
RadiusEdge=RadiusEdge(sortInd);
[V_innerEdge(:,1),V_innerEdge(:,2)]=pol2cart(ThetaEdge,RadiusEdge);

[Xr]=linspacen(V_innerEdge(:,1),V_outerEdge(:,1),nElements/2+1); Xr(end+1,:)=Xr(1,:);
[Yr]=linspacen(V_innerEdge(:,2),V_outerEdge(:,2),nElements/2+1); Yr(end+1,:)=Yr(1,:);

[Fs2,Vs2] = surf2patch(Xr,Yr,zeros(size(Xr)));
Vs2=Vs2(:,1:2);

V=[V_centralMesh;Vs2];
F=[F_centralMesh;Fs2+size(V_centralMesh,1)];
C=[ones(size(F_centralMesh,1),1); 2*ones(size(Fs2,1),1); ];

indEdge=((size(V,1)-size(Xr,1))+1):size(V,1);

%% Removing double points

[F,V,~,IND_IND]=mergeVertices(F,V);
indEdge=IND_IND(indEdge(1:end-1));

%Scaling radius
[ThetaMesh,RadiusMesh]=cart2pol(V(:,1),V(:,2));
RadiusMesh=r*RadiusMesh;
[V(:,1),V(:,2)]=pol2cart(ThetaMesh,RadiusMesh);

varargout{1}=F; 
varargout{2}=V; 
varargout{3}=C; 
varargout{4}=indEdge; 
 
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
