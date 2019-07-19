function [varargout]=regionTriMesh3D(varargin)

% function [F,V,boundaryInd]=regionTriMesh3D(regionCell,pointSpacing,resampleCurveOpt,interpMethod)
% ------------------------------------------------------------------------
% This function creates a 3D triangulation for the region specified in the
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
% See also: regionTriMesh2D
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2013/11/21: Created
% 2013/11/21: Added varargin input parsing
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        regionCell=varargin{1};
        pointSpacing=[];
        resampleCurveOpt=[];
        interpMethod=[];
    case 2
        regionCell=varargin{1};
        pointSpacing=varargin{2};
        resampleCurveOpt=[];
        interpMethod=[];
    case 3
        regionCell=varargin{1};
        pointSpacing=varargin{2};
        resampleCurveOpt=varargin{3};
        interpMethod=[];
    case 4
        regionCell=varargin{1};
        pointSpacing=varargin{2};
        resampleCurveOpt=varargin{3};
        interpMethod=varargin{4};
end

if isempty(pointSpacing)
    D=zeros(1,numel(regionCell));
    for q=1:1:numel(regionCell)
        D(q)=mean(diff(polyCurveLength(regionCell{1})));
    end
    pointSpacing=mean(D);
end

if isempty(resampleCurveOpt)
    resampleCurveOpt=0;
end
   
if isempty(interpMethod)
    interpMethod='natural';
end

%% Parse 3D regions and convert to planar 2D

% TO DO: get sizes and do it properly, hint cellfun(@(x) size(x,1),regionCell)
npSets=cellfun(@(x) size(x,1),regionCell);
npTotal=sum(npSets);

%Collect all points
Vt=zeros(npTotal,3);
startInd=1; 
endInd=npSets(1); 
for qCurve=1:1:numel(regionCell)
    Vs=regionCell{qCurve};    
    Vt(startInd:endInd,:)=Vs;    
    startInd=startInd+npSets(qCurve);
    if qCurve<numel(regionCell)
        endInd=startInd+npSets(qCurve+1)-1;
    end
end

Vt_mean=mean(Vt,1);
Vtc=Vt-Vt_mean(ones(size(Vt,1),1),:);

%Get orientation
[R_fit]=pointSetPrincipalDir(Vtc);

%Create planar 2D sets
Vtcr=Vt;
startInd=1; 
endInd=npSets(1); 
regionCell2D=regionCell;
for qCurve=1:1:numel(regionCell)    
    Vs=Vt(startInd:endInd,:);
    Vs=Vs-Vt_mean(ones(size(Vs,1),1),:);
    Vs=(R_fit'*Vs')';     
    regionCell2D{qCurve}=Vs(:,[1 2]);    
    Vtcr(startInd:endInd,:)=Vs;    
    startInd=startInd+npSets(qCurve);
    if qCurve<numel(regionCell)
        endInd=startInd+npSets(qCurve+1)-1;
    end
end

%% Use regionTriMesh2D to mesh region in 2D

[F,V,boundaryInd]=regionTriMesh2D(regionCell2D,pointSpacing,resampleCurveOpt,0);

%% Interpolate to obtain Z coordinates

interpFunc = scatteredInterpolant(Vtcr(:,[1 2]),Vtcr(:,3),interpMethod);
V(:,3)=interpFunc(V(:,[1 2]));
V=(R_fit*V')';
V=V+Vt_mean(ones(size(V,1),1),:);

%% Check normal directions against input curve(s)

% V_test=regionCell{1};
% V_test=V_test(1:2,:);
% size(V_test)
% size(V)
% [~,e]=minDist(V_test,V);
% E=patchEdges(F,1);
% if ~any(E(:,1)==e(1) & E(:,2)==e(2))
%    F=fliplr(F); %Invert if edge does not seem member of faces 
% end

%% Collect ouput
varargout{1}=F;
varargout{2}=V;
varargout{3}=boundaryInd;

%% Smoothen
% TR = triangulation(F,V);
% boundEdges = freeBoundary(TR);
% boundaryInd=unique(boundEdges(:));
% smoothPar.n=50;
% smoothPar.Tolerance=0.001;
% smoothPar.RigidConstraints=boundaryInd;
% smoothPar.smoothMethod='HC';
% [V]=tesSmooth(F,V,[],smoothPar);

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
