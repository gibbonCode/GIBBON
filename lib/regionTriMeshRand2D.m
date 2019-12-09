function [F,V,C,DT]=regionTriMeshRand2D(regionCell,pointSpacing,SD,resampleCurveOpt,plotOn)

% function [F,V,C,DT]=regionTriMesh2D(regionCell,pointSpacing,plotOn)
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
% kevinmoerman@hotmail.com
% 2013/14/08
%------------------------------------------------------------------------

%% CONTROL PARAMETERS
interpMethod='pchip';
closeLoopOpt=1;
minConnectivity=3; 

%% PLOT SETTINGS
if plotOn==1
    figColor='w'; figColorDef='white';
    fontSize=20;
    fAlpha=1;
    markerSize=20;
    faceColor='r';
end

%% Get region curves and define constraints for Delaunay tesselation

V=[]; %Curve points
C=[]; %Constraints
nss=0; %Number of points parameter for constraint index correction
for qCurve=1:1:numel(regionCell)
    %Get curve
    Vs=regionCell{qCurve};
    
    %Resample curve
    if resampleCurveOpt==1
        %Calculate required number of points for curve
        np=round(max(pathLength(Vs))./pointSpacing);
        [Vs]=evenlySampleCurve(Vs,np,interpMethod,closeLoopOpt);
    end
    
    %Collect curve points
    V=[V;Vs];
    
    %Create curve constrains
    ns= size(Vs,1);
    Cs=[(1:ns)' [2:ns 1]'];
    C=[C;Cs+nss];
    nss=nss+ns;
end

%% DEFINE INTERNAL MESH SEED POINTS

%Mesh point spacing for aproximately equilateral triangular mesh
pointSpacingXY=[pointSpacing pointSpacing.*0.5.*sqrt(3)];

%region extrema
maxVi=max(V(:,[1 2]),[],1);
minVi=min(V(:,[1 2]),[],1);
maxV=maxVi+pointSpacingXY+(4*SD);
minV=minVi-pointSpacingXY-(4*SD);
% maxV=maxVi+pointSpacingXY+rand(size(pointSpacingXY)).*pointSpacingXY;
% minV=minVi-pointSpacingXY-rand(size(pointSpacingXY)).*pointSpacingXY;

%Region "Field Of View" size
FOV=abs(maxV-minV);

%Calculate number of points in each direction
FOV_dev=FOV./pointSpacingXY;

npXY=round(FOV_dev);

%Get mesh of seed points
xRange=linspace(minV(1),maxV(1),npXY(1));
yRange=linspace(minV(2),maxV(2),npXY(2));
[X,Y]=meshgrid(xRange,yRange);

%Random offset
X=X+SD(1)*randn(size(X));
Y=Y+SD(2)*randn(size(Y));

%Offset mesh in X direction to obtain aproximate equilateral triangular mesh
X(1:2:end,:)=X(1:2:end,:)+(pointSpacing/2);

%Adding seed points to list
V_add=[X(:) Y(:)];
V=[V;V_add]; %Total point set of curves and seeds
numPointsIni=size(V,1);

%% DERIVE CONSTRAINED DELAUNAY TESSELATION

%Initial Delaunay triangulation
DT = delaunayTriangulation(V(:,1),V(:,2),C);
V=DT.Points;
F=DT.ConnectivityList;

%Remove poorly connected points associated with poor triangles
[~,IND_V]=patchIND(F,V);
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

if numPoints==numPointsPost
    warning('No points removed in contrained Delaunay triangulation. Possibly due to large pointSpacing with respect to curve size. Meshing skipped!');
    F=[];
    V=[];
else
    
    %% CONSTRAINED SMOOTHENING OF THE MESH
    
    TR = triangulation(F,V);
    boundEdges = freeBoundary(TR);
    boundaryInd=unique(boundEdges(:));
    
    smoothPar.LambdaSmooth=0.5;
    smoothPar.n=250;
    smoothPar.Tolerance=0.01;
    smoothPar.RigidConstraints=boundaryInd;
    [V]=tesSmooth(F,V,[],smoothPar);
    DT.Points=V;
    
    %% PLOTTING
    if plotOn==1
        cFigure;
        title('Smoothened triangulated mesh','FontSize',fontSize);
        xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize);zlabel('Z','FontSize',fontSize);
        hold on;
        patch('faces',F,'vertices',V,'FaceColor',faceColor,'FaceAlpha',fAlpha);
        plotV(V(boundaryInd',:),'b.','MarkerSize',markerSize);
        axis equal; view(2); axis tight;  set(gca,'FontSize',fontSize); grid on;
        drawnow;
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
