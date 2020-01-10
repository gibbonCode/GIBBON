function [varargout]=regionTriMesh2D(varargin)

% function [F,V]=regionTriMesh2D(regionCell,pointSpacing,resampleCurveOpt,plotOn)
% ------------------------------------------------------------------------
% This function creates a 2D triangulation for the region specified in the
% variable regionCell. The mesh aims to obtain a point spacing as defined
% by the input pointSpacing.
% The function output contains the triangular faces in F, the vertices in V
% and the per triangle region indices in regionInd. By setting plotOn to 0
% or 1 plotting can be switched on or off.
%
% More on the specification of the region:
% The input variable regionCell is a cell array containing all the boundary
% curves, e.g. for a two curve region 1 we would have something like
% regionSpec{1}={V1,V2} where V1 and V2 are the boundary curves. Multiple
% curves may be given here. The first curve should form the outer boundary
% of the entire region, the curves that follow should define holes inside
% this boundary and the space inside them is therefore not meshed.
%
% Kevin Mattheus Moerman
%
% 2013/14/08
% 2017/18/04 Added varargin style with defaults for missing parameters
% 2017/18/04 Fixed bug which caused accidental removal of interior points
% for multiregion meshes or meshes containing holes.
%------------------------------------------------------------------------

%%

switch nargin
    case 1
        regionCell=varargin{1};
        pointSpacing=[];
        resampleCurveOpt=0;
        plotOn=0;
        gridType='tri';
    case 2
        regionCell=varargin{1};
        pointSpacing=varargin{2};
        resampleCurveOpt=0;
        plotOn=0;
        gridType='tri';
    case 3
        regionCell=varargin{1};
        pointSpacing=varargin{2};
        resampleCurveOpt=varargin{3};
        plotOn=0;
        gridType='tri';
    case 4
        regionCell=varargin{1};
        pointSpacing=varargin{2};
        resampleCurveOpt=varargin{3};
        plotOn=varargin{4};
        gridType='tri';
    case 5
        regionCell=varargin{1};
        pointSpacing=varargin{2};
        resampleCurveOpt=varargin{3};
        plotOn=varargin{4};
        gridType=varargin{5};
end

if isempty(pointSpacing)
    %If pointSpacing is empty it will be based on the point spacing of the
    %input curves
    pointSpacingCurves=zeros(1,numel(regionCell));
    for q=1:1:numel(regionCell)
        V_now=regionCell{q};
        pointSpacingCurves(q)=mean(sqrt(sum(diff(V_now,1,1).^2,2)));        
    end    
    pointSpacing=mean(pointSpacingCurves);
end

%% CONTROL PARAMETERS
interpMethod='linear';
closeLoopOpt=1;
minConnectivity=4; %Minimum connectivity for points

%% PLOT SETTINGS
if plotOn==1
    fontSize=20;
    faceAlpha=1;
    markerSize=20;
    c=gjet(4);
    faceColor=c(4,:);
end

%% Get region curves and define constraints for Delaunay tesselation

V=[]; %Curve points
VSS=[]; %Curve points
C=[]; %Constraints
nss=0; %Number of points parameter for constraint index correction
for qCurve=1:1:numel(regionCell)
    %Get curve
    Vs=regionCell{qCurve};
    
    %Resample curve evenly based on point spacing
    np=ceil(max(pathLength(Vs))./pointSpacing); %Calculate required number of points for curve
    [Vss]=evenlySampleCurve(Vs,np,interpMethod,closeLoopOpt);
    
    %Create refined set for distance based edge point removal
    [Vss_split]=evenlySampleCurve(Vs,np*2,interpMethod,closeLoopOpt);
    
    %Resample curve
    if resampleCurveOpt==1
        Vs=Vss;
    end
    
    %Collect curve points
    V=[V;Vs]; %Original or interpolated set
    
    VSS=[VSS;Vss_split]; %Interpolated set
    
    %Create curve constrains
    ns= size(Vs,1);
    Cs=[(1:ns)' [2:ns 1]'];
    C=[C;Cs+nss];
    nss=nss+ns;
end

%% DEFINE INTERNAL MESH SEED POINTS

%region extrema
maxVi=max(V(:,[1 2]),[],1);
minVi=min(V(:,[1 2]),[],1);

switch gridType
    case 'tri'
        [~,Vt]=triMeshEquilateral(minVi,maxVi,pointSpacing);
        X=Vt(:,1); Y=Vt(:,2);
    case 'quad'
        w=maxVi-minVi;
        n=round(w./pointSpacing);
        [X,Y]=meshgrid(linspace(minVi(1),maxVi(1),n(1)),linspace(minVi(2),maxVi(2),n(2)));
        X=X(:); Y=Y(:);
end

%% Remove edge points

[D,~]=minDist([X Y],VSS);
L=D>(pointSpacing*0.5*sqrt(3))/2;

X=X(L);
Y=Y(L);

%%
%Adding seed points to list
V_add=[X(:) Y(:)];

V=[V;V_add]; %Total point set of curves and seeds
numPointsIni=size(V,1);

%% DERIVE CONSTRAINED DELAUNAY TESSELATION

%Initial Delaunay triangulation
DT = delaunayTriangulation(V(:,1),V(:,2),C);
V=DT.Points;
F=DT.ConnectivityList;

%% PLOTTING
    
%Remove poorly connected points associated with poor triangles
[~,IND_V]=tesIND(F,V,0); % [~,IND_V]=patchIND(F,V);

connectivityCount=sum(IND_V>0,2);
logicPoorConnectivity=connectivityCount<=minConnectivity;
logicConstraints=false(size(logicPoorConnectivity));
logicConstraints(C(:))=1;
logicRemoveList=logicPoorConnectivity & ~logicConstraints;
logicKeepList=~logicRemoveList;
V=V(logicKeepList,:); %Remove points
numPoints=size(V,1);

%Fix constraints list
indNew=(1:numPoints)';
indFix=zeros(numPointsIni,1);
indFix(logicKeepList)=indNew;
C=indFix(C);

%Redo Delaunary triangulation
DT = delaunayTriangulation(V(:,1),V(:,2),C);
V=DT.Points;
F=DT.ConnectivityList;

%Remove faces not inside region
L = isInterior(DT);
F=F(L,:);

%Removing unused points
indUni=unique(F(:));
indNew=nan(size(V,1),1);
indNew(indUni)=(1:numel(indUni))';
F=indNew(F);
V=V(indUni,:);

numPointsPost=size(V,1);

boundEdges = patchBoundary(F,V);
boundaryInd=unique(boundEdges(:));

if numPoints==numPointsPost
    warning('No points removed in contrained Delaunay triangulation. Possibly due to large pointSpacing with respect to curve size. Meshing skipped!');
%     F=[];
%     V=[];
else
    
    %% CONSTRAINED SMOOTHENING OF THE MESH
    
%     switch gridType
%         case 'tri'           
            smoothPar.LambdaSmooth=0.5;
            smoothPar.n=250;
            smoothPar.Tolerance=0.01;
            smoothPar.RigidConstraints=boundaryInd;
            [V]=tesSmooth(F,V,[],smoothPar);
%     end
    
    %% PLOTTING
    if plotOn==1
        cFigure;
        title('The meshed model','FontSize',fontSize);
        xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);zlabel('Z','FontSize',fontSize);
        hold on;
        gpatch(F,V,faceColor,'k',faceAlpha);
        plotV(V(boundaryInd',:),'b.','MarkerSize',markerSize);        
        axis equal; view(2); axis tight;  set(gca,'FontSize',fontSize); grid on;
        drawnow;
    end
    
end
 
%% Collect ouput
varargout{1}=F;
varargout{2}=V;
varargout{3}=boundaryInd;

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
