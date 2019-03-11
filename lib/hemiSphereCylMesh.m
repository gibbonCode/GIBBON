function [varargout]=hemiSphereCylMesh(inputStruct)

% function [Ft,Vt,Ct]=hemiSphereCylMesh(inputStruct)
% ------------------------------------------------------------------------
% Creates the faces (Ft) and vertices (Vt) for a triangular surface
% describing a hemisphere connected to a cylinder. The input is a structure
% array S containint the following control parameters:
%   S.sphereRadius => The radius of the hemi-spher portion
%   S.nRefineRegions => Number of |subtri| refinements for icosahedron
%   S.cylinderHeight => height of the cylinder part
%   S.cylinderStepSize => Aproximate node spacing for cylinder portion
%
% See also: |hemiSphereRegionMesh|
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2013/13/08
% ------------------------------------------------------------------------

%% Defining hemi-sphere portion
patchType=inputStruct.patchType;
switch patchType
    case {'tri','tri_slash'} %Triangular sphere mesh        
        [Fs,Vs]=hemiSphereMesh(inputStruct.nRefine,inputStruct.sphereRadius,0);        
        Ft=Fs;
    case 'quad' %Create quadrilateral sphere mesh        
        [Vs,Fs]=platonic_solid(2,inputStruct.sphereRadius);
        
        for q=1:1:(inputStruct.nRefine)
            [Fs,Vs]=subQuad(Fs,Vs,1);
            [T,P,R] = cart2sph(Vs(:,1),Vs(:,2),Vs(:,3));
            [Vs(:,1),Vs(:,2),Vs(:,3)]=sph2cart(T,P,ones(size(R)).*inputStruct.sphereRadius);
        end
        
        %Get top hemi-sphere
        LV=Vs(:,3)>(0-eps(0));
        LF=all(LV(Fs),2);        
        Fs=Fs(LF,:); %Crop faces
        
        %Remove unused vertices and fix face indices accordingly
        uniInd=unique(Fs(:));
        indAll=nan(1,size(Vs,1));
        indAll(uniInd)=1:numel(uniInd);
        Vs=Vs(uniInd,:);
        Fs=indAll(Fs);
        Ft=[Fs(:,1) Fs(:,3) Fs(:,2); Fs(:,1) Fs(:,4) Fs(:,3)];
        Ft=fliplr(Ft); %flip face orientation
end
Vs(:,3)=Vs(:,3)-max(Vs(:,3)); %Set max at zero
%Rotate upside down
R=euler2DCM([pi 0 0]);
Vs=(R*Vs')';

%% Defining cylindrical portion

%Cylinder parameters
cylinderHeight=inputStruct.cylinderHeight; %Cylinder height
cylinderStepSize=inputStruct.cylinderStepSize; %Edge length in Z direction for cylinder

% Find hemi-sphere edge
[Eb]=patchBoundary(Ft,Vs);
indCurve=edgeListToCurve(Eb);
indCurve=indCurve(1:end-1);
V_edge=Vs(indCurve,:);
if isPolyClockwise(V_edge)
   V_edge=flipud(V_edge);  
end
%Determine Z step sizes
if isempty(cylinderStepSize)
    cylinderStepSize=mean(sqrt(sum(diff(V_edge,1,1).^2,2))); %mean point spacing on edge
end
nStepsCylinder=round(cylinderHeight./cylinderStepSize);

if nStepsCylinder<2
    nStepsCylinder=2; 
end

if strcmp(patchType,'tri')
    nStepsCylinder=nStepsCylinder+iseven(nStepsCylinder);
end

controlParameterStruct.numSteps=nStepsCylinder;
controlParameterStruct.depth=cylinderHeight; 
controlParameterStruct.patchType=patchType; 
controlParameterStruct.dir=1;
controlParameterStruct.closeLoopOpt=1; 
[Fc,Vc]=polyExtrude(V_edge,controlParameterStruct);

%% MERGING MODEL PORTIONS

%Join element sets
[Ft,Vt,Ct]=joinElementSets({Fs,Fc},{Vs,Vc});

%Removing double vertices and fixing face indices
[Ft,Vt]=mergeVertices(Ft,Vt);

%% Collect output
varargout{1}=Ft;
varargout{2}=Vt;
if nargout==3
    varargout{3}=Ct;
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
