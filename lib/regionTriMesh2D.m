function [varargout]=regionTriMesh2D(varargin)

% function [F,V,boundaryInd]=regionTriMesh2D(inputStructure)
% function [F,V,boundaryInd]=regionTriMesh2D(regionCell,pointSpacing,resampleCurveOpt,plotOn);
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
%
% Kevin Mattheus Moerman
%
% 2013/14/08
% 2017/18/04 Added varargin style with defaults for missing parameters
% 2017/18/04 Fixed bug which caused accidental removal of interior points
% for multiregion meshes or meshes containing holes.
% 2020/05/05 Updated to handle single input structure as alternative
% 2020/05/05 Added interior "must points" option
% 2020/11/04 Improved handling of all must point boundary
%------------------------------------------------------------------------

%% Parse input

defaultInputStructure.regionCell=[];
defaultInputStructure.pointSpacing=[];
defaultInputStructure.resampleCurveOpt=0;
defaultInputStructure.plotOn=0;
defaultInputStructure.gridType='tri';
defaultInputStructure.mustPointsInner=[];
defaultInputStructure.mustPointsBoundary=[];
defaultInputStructure.smoothIterations=250;
defaultInputStructure.removeDistInterior=[];
defaultInputStructure.removeDistBoundary=[];
defaultInputStructure.SD=[];

switch nargin
    case 1
        singleInput=varargin{1};
        if isa(singleInput,'cell')
            inputStructure.regionCell=singleInput;
        elseif isa(singleInput,'struct')
            inputStructure=singleInput;
        else
            error(['For single input the first input should be either a cell or a structure, currenty input is of class: ',class(singleInput)]);
        end
    case 2
        inputStructure.regionCell=varargin{1};
        inputStructure.pointSpacing=varargin{2};
    case 3
        inputStructure.regionCell=varargin{1};
        inputStructure.pointSpacing=varargin{2};
        inputStructure.resampleCurveOpt=varargin{3};
    case 4
        inputStructure.regionCell=varargin{1};
        inputStructure.pointSpacing=varargin{2};
        inputStructure.resampleCurveOpt=varargin{3};
        inputStructure.plotOn=varargin{4};
    case 5
        inputStructure.regionCell=varargin{1};
        inputStructure.pointSpacing=varargin{2};
        inputStructure.resampleCurveOpt=varargin{3};
        inputStructure.plotOn=varargin{4};
        inputStructure.gridType=varargin{5};
    case 6
        inputStructure.regionCell=varargin{1};
        inputStructure.pointSpacing=varargin{2};
        inputStructure.resampleCurveOpt=varargin{3};
        inputStructure.plotOn=varargin{4};
        inputStructure.gridType=varargin{5};
        inputStructure.mustPointsInner=varargin{6};
    case 7
        inputStructure.regionCell=varargin{1};
        inputStructure.pointSpacing=varargin{2};
        inputStructure.resampleCurveOpt=varargin{3};
        inputStructure.plotOn=varargin{4};
        inputStructure.gridType=varargin{5};
        inputStructure.mustPointsInner=varargin{6};
        inputStructure.mustPointsBoundary=varargin{7};
    case 8
        inputStructure.regionCell=varargin{1};
        inputStructure.pointSpacing=varargin{2};
        inputStructure.resampleCurveOpt=varargin{3};
        inputStructure.plotOn=varargin{4};
        inputStructure.gridType=varargin{5};
        inputStructure.mustPointsInner=varargin{6};
        inputStructure.mustPointsBoundary=varargin{7};
        inputStructure.smoothIterations=varargin{8};
end

%Check optionStruct against default
[inputStructure]=structComplete(inputStructure,defaultInputStructure,1); %Complement provided with default if missing or empty

%Access structure components
regionCell=inputStructure.regionCell;
pointSpacing=inputStructure.pointSpacing;
resampleCurveOpt=inputStructure.resampleCurveOpt;
plotOn=inputStructure.plotOn;
gridType=inputStructure.gridType;
V_must_inner=inputStructure.mustPointsInner;
V_must_boundary=inputStructure.mustPointsBoundary;
smoothIterations=inputStructure.smoothIterations;
removeDistInterior=inputStructure.removeDistInterior;
removeDistBoundary=inputStructure.removeDistBoundary;
SD=inputStructure.SD(:)';

if isempty(regionCell)
    error('Empty regionCell provided');
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

if isempty(resampleCurveOpt)
    resampleCurveOpt=0;
end

if isempty(plotOn)
    plotOn=0;
end

if isempty(removeDistInterior)
    removeDistInterior=pointSpacing/2;
end

if isempty(removeDistBoundary)
    removeDistBoundary=(pointSpacing*0.5*sqrt(3))/2;
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
    if ~isempty(V_must_boundary)
        [~,indMust]=minDist(V_must_boundary,Vs);
        indMust=unique(indMust); %Force unique
    else
        indMust=[];
    end
    
    [Vss]=evenlySpaceCurve(Vs,pointSpacing,interpMethod,closeLoopOpt,indMust);

    %Create refined set for distance based edge point removal
    [Vss_split]=subCurve(Vss,1,1);
    
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
if isempty(SD)
    maxVi=max(V(:,[1 2]),[],1);
    minVi=min(V(:,[1 2]),[],1);
else
    maxVi=max(V(:,[1 2]),[],1)+(4*SD);
    minVi=min(V(:,[1 2]),[],1)-(4*SD);
end

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

%Random offset
if ~isempty(SD)
    X=X+SD(1)*randn(size(X));
    Y=Y+SD(2)*randn(size(Y));
end

%% Remove edge points

D=minDist([X Y],VSS);
L=D>removeDistBoundary;
X=X(L);
Y=Y(L);

if ~isempty(V_must_inner)
    D=minDist([X Y],V_must_inner);
    L=D>removeDistInterior;
    X=X(L);
    Y=Y(L);
end

%% Define additional point set

V_add=[X(:) Y(:)];

%%
%Adding seed points to list

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
% Add desired interior points
if ~isempty(V_must_inner)
    V=[V;V_must_inner];
    [V,~,ind2]=unique(V,'rows');
    C=ind2(C);
end

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

boundEdges = patchBoundary(F);
boundaryInd=unique(boundEdges(:));

indMustPointsInner=[];
if numPoints==numPointsPost
    warning('No points removed in contrained Delaunay triangulation. Possibly due to large pointSpacing with respect to curve size.');
    %     F=[];
    %     V=[];
else
    
    %% CONSTRAINED SMOOTHENING OF THE MESH
    
    if ~isempty(V_must_inner)
        [~,indMustPointsInner]=minDist(V_must_inner,V);
    end
    
    if smoothIterations>0
        smoothPar.Method='LAP';
        smoothPar.n=smoothIterations;
        smoothPar.Tolerance=0.01;
        smoothPar.RigidConstraints=[boundaryInd(:); indMustPointsInner];        
        [V]=tesSmooth(F,V,[],smoothPar);
    end
    
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
if nargout>3
    varargout{4}=indMustPointsInner;
end
if nargout>4
    if ~isempty(V_must_boundary)
        [~,indMin]=minDist(V_must_boundary,V(boundaryInd,:));
        varargout{5}=boundaryInd(indMin);
    else
        varargout{5}=[];
    end
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
