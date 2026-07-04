function [Fn,Vn]=patchExtend(varargin)

% function [Fn,Vn]=patchExtend(F,V,Eb,optionStruct)
% ------------------------------------------------------------------------
% This function extends the input surface geometry allong the edges Eb
% using the local (or averaged) surface direction at these edges. 
%
% The input consists of: the input surface (patch) data, defined by the
% faces F and vertices V, the boundary edges Eb to extend the surface
% over, and an optionStruct structure containing all remaining and optional
% inputs. 
%
% The input structure may contain the following fields: 
%
% optionStruct.extendDistance -> %The distance to extend the surface by
%
% optionStruct.numSteps -> The number of steps to use in the extension
% direction, default is the distance divided by the average edge length. 
%
% optionStruct.plotOn -> option to turn on/off plotting, default is 0
%
% optionStruct.meshType -> the mesh type ('tri' or 'quad'), default is
% 'tri'. 
%
% optionStruct.globalDirection -> the global (should be seen as mean)
% direction for extension, default is empty ([]) and is instead based on
% the input mesh directions
%
% optionStruct.extendMethod=1; %Method to use for extruding (see below),
% default is 1. 
%
% Five extend methods have been implemented (default is 1): 
%
% extendMethod=1 -> Equal offset allong local direction
% extendMethod=2 -> Equal offset wrst mean direction allong local direction
% extendMethod=3 -> Varying offset allong local ending planar wrt mean direction
% extendMethod=4 -> Equal offset allong mean direction
% extendMethod=5 -> Equal offset allong mean direction ending planar wrt mean direction 
% 
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% Created: 2020/04/01
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        Eb=[];
        optionStruct=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        Eb=varargin{3};
        optionStruct=[];
    case 4
        F=varargin{1};
        V=varargin{2};
        Eb=varargin{3};
        optionStruct=varargin{4};
end

%Check for empty boundary edges
if isempty(Eb)
    Eb=patchBoundary(F,V); %use all boundary edges
end

%Check boundary point spacing
pointSpacing=mean(patchEdgeLengths(Eb,V));

%Create default option structure
defaultOptionStruct.numSteps=[];
defaultOptionStruct.plotOn=0;
defaultOptionStruct.meshType=[];
defaultOptionStruct.globalDirection=[];
defaultOptionStruct.extendDistance=[]; %Extend distance
defaultOptionStruct.extendMethod=1; %Method to use for extruding

%Parse option structure
if isempty(optionStruct)
    optionStruct=defaultOptionStruct;
else
    %Check optionStruct against default
    [optionStruct]=structComplete(optionStruct,defaultOptionStruct,1); %Complement provided with default if missing or empty
end

%Access input parameters
numSteps=optionStruct.numSteps;
plotOn=optionStruct.plotOn;
meshType=optionStruct.meshType;
ns_mean=optionStruct.globalDirection; 
extendDistance=optionStruct.extendDistance;
extendMethod=optionStruct.extendMethod;

if isempty(extendDistance)
    extendDistance=pointSpacing;
end

if isempty(numSteps)
    numSteps=ceil(extendDistance./pointSpacing);    
end
if numSteps<2
    numSteps=2; 
end

if isempty(meshType)
    switch size(F,2)
        case 3
            meshType='tri';            
        case 4
            meshType='quad';
        otherwise
            meshType='quad';
    end
end

%%

indEdgePoints=unique(Eb(:)); %Vertex indices of all edge points

[~,~,N]=patchNormal(F,V); %Get surface vertex normal vectors
n=N((indEdgePoints),:); %Surface normal directions at edge start points
[Eb,V1]=patchCleanUnused(Eb,V);
numEdgePoints=size(V1,1);

[~,~,ne]=edgeVec(Eb,V1); %Get edge vectors
ne=vecnormalize(ne); %Edge vectors

ns=vecnormalize(cross(ne,n)); %Allong surface extrude directions orthogonal to both
if isempty(ns_mean)
    ns_mean=vecnormalize(mean(ns,1)); %Mean extrude direction
end

switch extendMethod
    case 1 %Equal offset allong local direction
        V2=V1+extendDistance*ns;
    case 2 %Equal offset wrst mean direction allong local direction
        f=extendDistance./dot(ns,ns_mean(ones(size(ns,1),1),:),2);
        V2=V1+f.*ns;
    case 3 %Varying offset allong local ending planar wrt mean direction
        d=dot(ns,ns_mean(ones(size(ns,1),1),:),2);
        f=extendDistance./d;
        dd=dot(V1-mean(V1,1),ns_mean(ones(numEdgePoints,1),:),2);
        [~,indMax]=max(dd);
        dp=dot(V1-V1(indMax,:),ns_mean(ones(numEdgePoints,1),:),2);
        fp=dp./d;
        V2=V1-fp.*ns+f.*ns;
    case 4 %Equal offset allong mean direction
        V2=V1+extendDistance*ns_mean(ones(numEdgePoints,1),:);
    case 5 %Equal offset allong mean direction ending planar wrt mean direction   
        dd=dot(V1-mean(V1,1),ns_mean(ones(numEdgePoints,1),:),2);
        [~,indMax]=max(dd);
        fp=dot(V1-V1(indMax,:),ns_mean(ones(numEdgePoints,1),:),2);
        V2=V1-fp.*ns_mean(ones(numEdgePoints,1),:)+extendDistance*ns_mean(ones(numEdgePoints,1),:);
    otherwise
        error('Invalid extrude option provided, use 1-5');
end

%Create quadrilateral faces
X=linspacen(V1(:,1),V2(:,1),numSteps);
Y=linspacen(V1(:,2),V2(:,2),numSteps);
Z=linspacen(V1(:,3),V2(:,3),numSteps);

Vn=[X(:) Y(:) Z(:)]; %Vertex set
Fn=repmat([fliplr(Eb) Eb],numSteps-1,1);
q=repmat(0:(numSteps-2),size(Eb,1),1);
q=q(:);
Q=[q q q+1 q+1]*size(V1,1);
Fn=Fn+Q;

switch meshType
    case 'tri'
        [Fn,Vn]=quad2tri(Fn,Vn,'a');
    case 'quad'
        
end

%%
if plotOn==1
    cFigure; 
    subplot(1,2,1); hold on;
    title(['Extend method: ',num2str(extendMethod)]);
    hp1(1)=gpatch(F,V,'kw','none',0.5);
    hp1(2)=gpatch(Fn,Vn,'bw','b',0.5,1);    
    hp1(3)=quiverVec(V1,n,pointSpacing,'r');
    hp1(4)=quiverVec(V1,ne,pointSpacing,'g');
    hp1(5)=quiverVec(V1,ns,pointSpacing,'b');    
    hp1(6)=quiverVec(mean(V1,1),ns_mean,extendDistance+5*pointSpacing,'k');
    legend(hp1,{'Input surface','Extended surface','Vertex surface normal',...
               'vertex edge vector','vertex extrude vector','Global or mean direction'},...
               'Location','SouthOutSide');    
    axisGeom; camlight headlight;   
        
    subplot(1,2,2); hold on;    
    hp2(1)=gpatch(F,V,'kw','k',1);
    hp2(2)=gpatch(Fn,Vn,'bw','k',1);        
    legend(hp2,{'Input surface','Extended surface'},'Location','SouthOutSide');        
    axisGeom; camlight headlight;   
    
    gdrawnow;    
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
