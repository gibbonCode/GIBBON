function [F,V]=polyExtrude(Vc,cPar)

%% COMPUTE CURVE LENGTH
D=max(pathLength(Vc)); %Compute curve length for point sampling
numPoints=size(Vc,1);

%% Cope with 2D input
if size(Vc,2)==2
    Vc(:,3)=0;
end

%% PARSE cPar

%Check what mode
if isfield(cPar,'pointSpacing') && isfield(cPar,'numSteps')
    error('Either specify pointSpacing or numSteps, not both.');
end

%Check extrudeMode, numSteps and pointSpacing settings
if isfield(cPar,'numSteps')
    extrudeMode=1;
else
    extrudeMode=0;
    if ~isfield(cPar,'pointSpacing')
        cPar.pointSpacing=[]; %Default is average point spacing
    end
    if isempty(cPar.pointSpacing)
        cPar.pointSpacing=D/numPoints; %Default is average point spacing
    end
end

%Check depth
if ~isfield(cPar,'depth')
    error('cPar.depth was not specified.');
end

%Check dir
if ~isfield(cPar,'dir')
    cPar.dir=0; %Default symmetric
end

%Check patchType
if ~isfield(cPar,'patchType')
    cPar.patchType='tri'; %Default triangles
end

%Check direction vector
if ~isfield(cPar,'n')    
    [R_fit]=pointSetPrincipalDir(Vc);
    nDir=R_fit(:,3);
    if dot(nDir,[0 0 1])<1 %Make the z direction the default upward direction
        nDir=-nDir;
    end
    cPar.n=nDir; %Default normal direction to polygon
end
cPar.n=vecnormalize(cPar.n);
cPar.n=cPar.n(:)';

%Check closeLoopOpt
if ~isfield(cPar,'closeLoopOpt')
    cPar.closeLoopOpt=0; %Default off
end

%% Extruding the skethed profile

switch extrudeMode
    case 1 %Use number of points
        
    case 0 %Resample curve and use pointSpacing
        %Set point spacing
        cPar.numSteps=round(cPar.depth./cPar.pointSpacing);
        
        %Resampling sketch with desired spacing
        D=max(pathLength(Vc)); %Computer curve length for point sampling
        n=round(D./cPar.pointSpacing); %Determine number of points based on curve length
        interpMethod='linear';
        [Vc]=evenlySampleCurve(Vc,n,interpMethod,cPar.closeLoopOpt); %Resampling curve
end

%Create coordinates in extrusion direction
switch cPar.dir
    case 0
        V_add=(cPar.depth/2).*cPar.n;
        Vc_start=Vc-V_add(ones(size(Vc,1),1),:);
        Vc_end=Vc+V_add(ones(size(Vc,1),1),:);
    case 1
        Vc_start=Vc;
        V_add=cPar.depth.*cPar.n;
        Vc_end=Vc+V_add(ones(size(Vc,1),1),:);
    case -1
        Vc_start=Vc;
        V_add=cPar.depth.*cPar.n;
        Vc_end=Vc-V_add(ones(size(Vc,1),1),:);
end

[F,V]=polyLoftLinear(Vc_start,Vc_end,cPar);
 
%% <-- GIBBON footer text --> 
% 
%     GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
%     image segmentation, image-based modeling, meshing, and finite element
%     analysis.
% 
%     Copyright (C) 2017  Kevin Mattheus Moerman
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
