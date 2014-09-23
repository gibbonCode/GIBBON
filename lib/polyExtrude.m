function [F,V]=polyExtrude(Vc,cPar)

%% COMPUTER CURVE LENGTH
D=max(pathLength(Vc)); %Computer curve length for point sampling
numPoints=size(Vc,1); 

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
if ~isfield(cPar,'depth');
    error('cPar.depth was not specified.');
end

%Check extrusion mode
if ~isfield(cPar,'dir');
    cPar.dir=0; %Default symmetric
end

%Check extrusion mode
if ~isfield(cPar,'patchType');
    cPar.patchType='tri'; %Default triangles
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
        interpMethod='pchip';
        closeLoopOpt=0;
        [Vc]=evenlySampleCurve(Vc,n,interpMethod,closeLoopOpt); %Resampling curve
end

%Create z coordinates in extrusion direction
switch cPar.dir
    case 0
        z=linspace(-cPar.depth/2,cPar.depth/2,cPar.numSteps);
    case 1
        z=linspace(0,cPar.depth,cPar.numSteps);
    case -1
        z=linspace(-cPar.depth,0,cPar.numSteps);        
end

%Create X Y Z quad meshes
X=repmat(Vc(:,1),[1,cPar.numSteps]);
Z=repmat(Vc(:,2),[1,cPar.numSteps]);
Y=repmat(z(:),[1,size(Vc,1)])';

%Create quad patch data
[F,V,~] = surf2patch(X,Y,Z);

switch cPar.patchType
    case 'tri' %Convert quads to triangles        
        F=[F(:,1) F(:,3) F(:,2); F(:,1) F(:,4) F(:,3)];
end

