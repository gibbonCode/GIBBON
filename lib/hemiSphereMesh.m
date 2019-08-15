function [varargout]=hemiSphereMesh(varargin)

% function [F,V,C]=hemiSphereMesh(nRefineSteps,sphereRadius,closeOpt)
%-------------------------------------------------------------------------
% Creates the patch data (faces F, vertices V, and colordata C) for a
% hemisphere with a radius equal to sphereRadius. The hemisphere is
% refined nRefineSteps times and is closed if closeOpt==1
%
%
% 
%-------------------------------------------------------------------------

%% Parse input

switch nargin
    case 0
        nRefineSteps=0;
        sphereRadius=1;
        closeOpt=0;        
    case 1
        nRefineSteps=varargin{1};
        sphereRadius=1;
        closeOpt=0;
    case 2
        nRefineSteps=varargin{1};
        sphereRadius=varargin{2};
        closeOpt=0;
    case 3
        nRefineSteps=varargin{1};
        sphereRadius=varargin{2};
        closeOpt=varargin{3};
end

%% Create geodesic dome

%Create a geodesic dome by sub-triangulating an icosahedron
[V,F]=platonic_solid(4,sphereRadius); 
for q=1:1:(nRefineSteps+1)
    [F,V]=subtri(F,V,1);
    [az,elev,r] =cart2sph(V(:,1),V(:,2),V(:,3));
    [V(:,1),V(:,2),V(:,3)] =sph2cart(az,elev,sphereRadius.*ones(size(r)));    
end

%Rotate icosahedron so cut is clean
v3=vecnormalize(V(3,:)-V(1,:));
v1=[1 0 0];
v2=cross(v3,v1);
R=[v1' v2' v3'];
V=(R*V')';

%% Get half of sphere
[VF]=patchCentre(F,V); %Get face centre coordinates
logicKeep=VF(:,3)>0; %Logic to select faces above zero
F=F(logicKeep,:); %The cropped face set
[F,V]=patchCleanUnused(F,V); %Remove unused vertices
C=(1:1:size(F,1))'; %Color data initialized as face indices

%% Close bottom if needed
if closeOpt==1
   [Eb]=patchBoundary(F,V);
   indCurve=edgeListToCurve(Eb); 
   indCurve=indCurve(1:end-1);
   [Fb,Vb]=regionTriMesh2D({V(indCurve,[1 2])},[],0,0); %2D meshing of disc
   Vb(:,3)=mean(V(indCurve,3)); %Set z-coordinate
   Fb=fliplr(Fb); %Point normal downward
   [F,V,C]=joinElementSets({F,Fb},{V,Vb}); %Join element sets
   [F,V]=mergeVertices(F,V); %Merge nodes
end

%% Collect output

varargout{1}=F;
varargout{2}=V;
if nargout>2
    varargout{3}=C;
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
