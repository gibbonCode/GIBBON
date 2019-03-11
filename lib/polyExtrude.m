function [varargout]=polyExtrude(Vc,controlParameterStruct)

%% COMPUTE CURVE LENGTH
D=max(pathLength(Vc)); %Compute curve length for point sampling
numPoints=size(Vc,1);

%% Cope with 2D input
if size(Vc,2)==2
    Vc(:,3)=0;
end

%% Parse input

%Check what mode
if isfield(controlParameterStruct,'pointSpacing') && isfield(controlParameterStruct,'numSteps')
    error('Either specify pointSpacing or numSteps, not both.');
end

%Check extrudeMode, numSteps and pointSpacing settings
if isfield(controlParameterStruct,'numSteps')
    extrudeMode=1;
else
    extrudeMode=0;
    if ~isfield(controlParameterStruct,'pointSpacing')
        controlParameterStruct.pointSpacing=[]; %Default is average point spacing
    end
    if isempty(controlParameterStruct.pointSpacing)
        controlParameterStruct.pointSpacing=D/numPoints; %Default is average point spacing
    end
end

%Check depth
if ~isfield(controlParameterStruct,'depth')
    error('cPar.depth was not specified.');
end

%Check dir
if ~isfield(controlParameterStruct,'dir')
    controlParameterStruct.dir=0; %Default symmetric
end

%Check patchType
if ~isfield(controlParameterStruct,'patchType')
    controlParameterStruct.patchType='tri'; %Default triangles
end

%Check direction vector
if ~isfield(controlParameterStruct,'n')    
    [R_fit]=pointSetPrincipalDir(Vc);
    nDir=R_fit(:,3);
    if dot(nDir,[0 0 1])<0 %Make the z direction the default upward direction
        nDir=-nDir;        
    end
    controlParameterStruct.n=nDir; %Default normal direction to polygon
end
controlParameterStruct.n=vecnormalize(controlParameterStruct.n);
controlParameterStruct.n=controlParameterStruct.n(:)';

%Check closeLoopOpt
if ~isfield(controlParameterStruct,'closeLoopOpt')
    controlParameterStruct.closeLoopOpt=0; %Default off
end

%% Prepare for extrude

switch extrudeMode
    case 0 %Resample curve and use pointSpacing
        %Set point spacing
        controlParameterStruct.numSteps=round(controlParameterStruct.depth./controlParameterStruct.pointSpacing);
        
        %Resampling sketch with desired spacing
        D=max(pathLength(Vc)); %Computer curve length for point sampling
        n=round(D./controlParameterStruct.pointSpacing); %Determine number of points based on curve length
        interpMethod='linear';
        [Vc]=evenlySampleCurve(Vc,n,interpMethod,controlParameterStruct.closeLoopOpt); %Resampling curve
    case 1 %Use number of points
        
end

%Create coordinates in extrusion direction
switch controlParameterStruct.dir
    case 0
        V_add=(controlParameterStruct.depth/2).*controlParameterStruct.n;
        Vc_start=Vc-V_add(ones(size(Vc,1),1),:);
        Vc_end=Vc+V_add(ones(size(Vc,1),1),:);
    case 1
        Vc_start=Vc;
        V_add=controlParameterStruct.depth.*controlParameterStruct.n;
        Vc_end=Vc+V_add(ones(size(Vc,1),1),:);
    case -1
        Vc_start=Vc;
        V_add=controlParameterStruct.depth.*controlParameterStruct.n;
        Vc_end=Vc-V_add(ones(size(Vc,1),1),:);
end

%% Extruding the skethed profile using polyLoftLinear

[F,V,indStart,indEnd]=polyLoftLinear(Vc_start,Vc_end,controlParameterStruct);
varargout{1}=F;
varargout{2}=V;
varargout{3}=indStart;
varargout{4}=indEnd;

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
