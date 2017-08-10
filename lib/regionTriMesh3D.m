function [F,V]=regionTriMesh3D(regionCell,pointSpacing,resampleCurveOpt,interpMethod)

% function [F,V]=regionTriMesh3D(regionCell,pointSpacing,resampleCurveOpt,interpMethod)
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
% kevinmoerman@hotmail.com
% 2013/21/11
%------------------------------------------------------------------------

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

[F,V]=regionTriMesh2D(regionCell2D,pointSpacing,resampleCurveOpt,0);

%% Interpolate to obtain Z coordinates

interpFunc = scatteredInterpolant(Vtcr(:,[1 2]),Vtcr(:,3),interpMethod);
V(:,3)=interpFunc(V(:,[1 2]));
V=(R_fit*V')';
V=V+Vt_mean(ones(size(V,1),1),:);

%%

TR = triangulation(F,V);
boundEdges = freeBoundary(TR);
boundaryInd=unique(boundEdges(:));

smoothPar.n=50;
smoothPar.Tolerance=0.001;
smoothPar.RigidConstraints=boundaryInd;
smoothPar.smoothMethod='HC';
[V]=tesSmooth(F,V,[],smoothPar);



 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
