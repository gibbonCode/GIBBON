function [VV]=triPatchSmoothRayTraced(F,V,optionStruct)

nMax=optionStruct.n;
LambdaSmooth=0.5;
indRigid=optionStruct.RigidConstraints;

C=patchConnectivity(F,V,{'vv','vf'});

IND_V=C.vertex.vertex;
IND_F=C.vertex.face;

logicValid=IND_V>0;

logicValidFaces=IND_F>0;
[I,~]=ind2sub(size(IND_F),find(logicValidFaces));
indFaces=IND_F(logicValidFaces);

%Ray trace options
rayTraceStruct.eps      = 1e-6;
rayTraceStruct.triangle = 'one sided';
rayTraceStruct.ray      = 'ray';
rayTraceStruct.border   = 'normal';

nDims=size(V,2); %Number of dimensions

PP=V;
VV=V;

for qIter=1:nMax
    P=VV;
    
    %% Laplacian smoothing step
    % Calculate Laplacian smoothed set
    for qDim=1:1:nDims %Loop for all dimensions
        Xp=NaN(size(IND_V));
        Xp(logicValid)=P(IND_V(logicValid),qDim);
        Xp=mean(Xp,2,'omitnan');
        PP(:,qDim)=Xp;
    end
    P=P+LambdaSmooth.*(PP-P);
    
    %% Ray tracin step
    % Push Laplacian smoothed coordinates back to previous surface state
    
    % Get current normals
    [~,~,N]=patchNormal(F,P); 
    Vi = triangleRayIntersection(P(I,:),-N(I,:),VV,F(indFaces,:),rayTraceStruct);
    
    for qDim=1:1:nDims %Loop for all dimensions
        Xp=NaN(size(IND_F));
        Xp(logicValidFaces)=Vi(:,qDim);
        Xp=mean(Xp,2,'omitnan'); %Use mean if multiple intersections found
        Xp(isnan(Xp))=VV(isnan(Xp),qDim);
        P(:,qDim)=Xp;
    end
   
    %% Put back constrained points
    if ~isempty(indRigid)
        P(indRigid,:)=V(indRigid,:);
    end
    
    %% Fix nans
    logicNan=any(isnan(P),2);
    P(logicNan,:)=VV(logicNan,:);
    VV=P;
    
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
