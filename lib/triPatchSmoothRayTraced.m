function [VV]=triPatchSmoothRayTraced(F,V,optionStruct)

nMax=optionStruct.n;
LambdaSmooth=0.5;
indRigid=[];

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


