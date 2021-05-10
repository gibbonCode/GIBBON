function [FLs,VLs,logicAttachFaces,indCurve,Eb]=localRefineMap(FL,VL,V_SLIL_Lunate,distGrow,numGrowStepsRefine,cParSmoothCurve)

% Refine locally

[~,indClose]=minDist(V_SLIL_Lunate,VL);
logicAttachVertices=false(size(VL,1),1);
logicAttachVertices(indClose)=1;
logicAttachFaces=any(logicAttachVertices(FL),2);
logicAttachVertices=false(size(VL,1),1);

for q=1:1:numGrowStepsRefine
    logicAttachVertices(FL(logicAttachFaces))=1;
    logicAttachFaces=any(logicAttachVertices(FL),2);
end

% inputStruct.F=FL;
% inputStruct.V=VL;
% inputStruct.indFaces=find(logicAttachFaces);
% [outputStruct]=subTriLocal(inputStruct);
% FLs=outputStruct.F;
% VLs=outputStruct.V;

[FLs,VLs]=subTriDual(FL,VL,logicAttachFaces);

%%
% Map attachement site
[~,indClose]=minDist(V_SLIL_Lunate,VLs);
logicAttachVertices=false(size(VLs,1),1);
logicAttachVertices(indClose)=1;

[DLs]=minDist(VLs,VLs(logicAttachVertices,:));
logicAttachVertices=DLs<=distGrow;

logicAttachFaces=any(logicAttachVertices(FLs),2);
[logicAttachFaces]=triSurfLogicSharpFix(FLs,logicAttachFaces,3);

Eb=patchBoundary(FLs(logicAttachFaces,:),VLs);
indCurve=edgeListToCurve(Eb);
indCurve=indCurve(1:end-1);

logicAttachFaces_temp=logicAttachFaces;
indAttachVertices=unique(FLs(logicAttachFaces_temp,:));

logicHoldOn=true(size(VLs,1),1);
logicHoldOn(indAttachVertices)=0;
cParSmoothCurve.RigidConstraints=find(logicHoldOn);

[VLs]=patchSmooth(FLs(logicAttachFaces,:),VLs,[],cParSmoothCurve);

end

%%
function [F_ligament,V_ligament]=fixOverlap(F_ligament,V_ligament,VL,distanceTooClose,stepSizeTooCloseFix)

indBoundary=unique(patchBoundary(F_ligament,V_ligament));

D=minDist(V_ligament,VL);
logicClose=D<distanceTooClose;
logicClose(indBoundary)=0;
if any(logicClose)
    while 1
        [~,~,Nv]=patchNormal(F_ligament,V_ligament);
        D=minDist(V_ligament,VL);
        logicClose=D<distanceTooClose;
        logicClose(indBoundary)=0;
        if nnz(logicClose)==0
            break
        end
        V_ligament(logicClose,:)=V_ligament(logicClose,:)+stepSizeTooCloseFix*Nv(logicClose,:);
    end
    
    %Smooth offset result
    controlParameterSmooth.RigidConstraints=indBoundary;
    controlParameterSmooth.Method='HC';
    controlParameterSmooth.n=15;
    [V_ligament]=patchSmooth(F_ligament,V_ligament,[],controlParameterSmooth);
end

end
