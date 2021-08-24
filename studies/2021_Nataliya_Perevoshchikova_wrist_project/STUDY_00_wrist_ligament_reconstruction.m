%% STUDY_00_wrist_ligament_reconstruction
% Below is a demonstration for:
%
% * Reconstruction of scapholunate dorsal ligament in the wrist from
% attachments sites on scaphoid and lunate bones
% * Create solid mesh using tetgen
% * Importing and visualizing results
% * Filling up with fibers alone certain direction

%% Keywords
% * ligament, tendon
% * meshing
% * fibers


clear; close all; clc;

% Plot settings
fontSize=30;

%% Control parameters
%Parameters for SLIL reconstruction
numSmoothStepsMain=5;
numSmoothStepsCurve=15;
numSmoothStepsLoft=50;
numSmoothStepsTransition=10;
nRefineAttachmentSite=4;
numGrowStepsRefine=1;
csapsSmoothPar=0.99; %Cubic smoothening spline smoothening parameter
distGrow=0.21;
stepSizeTooCloseFix=0.1;
distanceTooClose=0.1;
fixOverlapOpt=0;
viewSpec = [-4 15 -5];

cParSmoothMain.n=numSmoothStepsMain;
cParSmoothMain.Method='HC';

cParSmoothCurve.n=numSmoothStepsCurve;
cParSmoothCurve.Method='HC';

cParSmoothLoft.n=numSmoothStepsLoft;
cParSmoothLoft.Method='LAP';

cParSmoothTransition.n=numSmoothStepsTransition;
cParSmoothTransition.Method='HC';


% Surface model file names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
loadPathSurfaces=fullfile(defaultFolder,'2021_Nataliya_Perevoshchikova_wrist_project','data','STL');

%Load geometries
matname = fullfile(loadPathSurfaces,'Temp.mat');
load(matname);


%Visualization
hf=cFigure;
hold on;

hl(1)=gpatch(FL,VL,'g','none',0.5);
hl(2)=gpatch(FS,VS,'y','none',0.5);

hl(3)=gpatch(F_SLIL_Lunate,V_SLIL_Lunate,'rw','k',1);
hl(4)=gpatch(F_SLIL_Scaphoid,V_SLIL_Scaphoid,'rw','k',1);

legend(hl,{'Lunate','Scaphoid','Ligament attachment sites'},'Position',[0.2 0.52 0.15 0.0869]);
grid on; axis equal;
axisGeom(gca,fontSize);
camlight headlight;
view(viewSpec);
drawnow; 

%% Local refinement and mapping of attachement site
[FL,VL,CL,indCurve_L,Eb_L]=localRefineMap(FL,VL,V_SLIL_Lunate,distGrow,numGrowStepsRefine,cParSmoothCurve);
[FS,VS,CS,indCurve_S,Eb_S]=localRefineMap(FS,VS,V_SLIL_Scaphoid,distGrow,numGrowStepsRefine,cParSmoothCurve);

%%
nL=vecnormalize(mean(patchNormal(FL(CL==1,:),VL),1));
nS=vecnormalize(mean(patchNormal(FS(CS==1,:),VS),1));

%%
% Visualize raw attachement sites

cFigure; hold on;
plotV(VL(indCurve_L,:),'g.-','MarkerSize',25,'LineWidth',3);
plotV(VS(indCurve_S,:),'g.-','MarkerSize',25,'LineWidth',3);

gpatch(FL,VL,'w','none',0.5);
gpatch(FS,VS,'w','none',0.5);

gpatch(F_SLIL_Lunate,V_SLIL_Lunate,'r');
gpatch(F_SLIL_Scaphoid,V_SLIL_Scaphoid,'r');
camlight headlight;
axisGeom(gca,fontSize);
view(viewSpec);
drawnow;

%% Edge splitting on the meshes to create same number of points on each boundary curve
while 1
    logicBoundary_L=false(size(VL,1),1);
    logicBoundary_L(indCurve_L)=1;
    
    logicBoundary_S=false(size(VS,1),1);
    logicBoundary_S(indCurve_S)=1;
    
    numPoints_S=numel(indCurve_S);
    numPoints_L=numel(indCurve_L);
    numSplit=abs(numPoints_S-numPoints_L); %Number of split iterations
    
    if numSplit==0
        break
    else
        if numPoints_L<numPoints_S
            D=sqrt(sum(diff(VL(indCurve_L,:),1,1).^2,2));
            [~,indMax]=max(D);
            edgeNow=[indCurve_L(indMax) indCurve_L(indMax+1)];
            
%             [Fs,Vs,Cs,CVs,CFs]=triEdgeSplit(F,V,E,CF,CV)
            [FL,VL,~,CL,CV]=triEdgeSplit(FL,VL,edgeNow,CL,logicBoundary_L);
            
            indBoundary=find(CV>0);
            
            E=patchEdges(FL,1);
            EL=E(all(ismember(E,indBoundary),2),:);
            
            [indCurve_L]=edgeListToCurve(EL);
            indCurve_L=indCurve_L(1:end-1);
            
        elseif numPoints_S<numPoints_L
            D=sqrt(sum(diff(VS(indCurve_S,:),1,1).^2,2));
            [~,indMax]=max(D);
            edgeNow=[indCurve_S(indMax) indCurve_S(indMax+1)];
            
            [FS,VS,~,CS,CV]=triEdgeSplit(FS,VS,edgeNow,CS,logicBoundary_S);
            indBoundary=find(CV>0);
            
            E=patchEdges(FS,1);
            E_tendon_contour_cut=E(all(ismember(E,indBoundary),2),:);
            
            [indCurve_S]=edgeListToCurve(E_tendon_contour_cut);
            indCurve_S=indCurve_S(1:end-1);
        end
    end
end

%%
% Visualization
cFigure; hold on;
plotV(VL(indCurve_L,:),'g.-','MarkerSize',25,'LineWidth',3);
plotV(VS(indCurve_S,:),'g.-','MarkerSize',25,'LineWidth',3);

gpatch(FL,VL,'w','none',0.5);
gpatch(FS,VS,'w','none',0.5);

gpatch(F_SLIL_Lunate,V_SLIL_Lunate,'r');
gpatch(F_SLIL_Scaphoid,V_SLIL_Scaphoid,'r');
camlight headlight;
axisGeom(gca,fontSize);
view(viewSpec);
drawnow;

%% Reorder curves
[~,indMin]=minDist(VL(indCurve_L(1),:),VS(indCurve_S,:));

if indMin>1
    indCurve_S=[indCurve_S(indMin:end) indCurve_S(1:indMin-1)];
end

D1=sum(sqrt(sum((VL(indCurve_L,:)-VS(indCurve_S,:)).^2,2)));
D2=sum(sqrt(sum((VL(indCurve_L,:)-VS(flip(indCurve_S),:)).^2,2)));

if D2<D1
    indCurve_S=flip(indCurve_S);
end

%% Construction of Ligament 
[~,~,NLv]=patchNormal(FL,VL);
[~,~,NSv]=patchNormal(FS,VS);

ds=mean(patchEdgeLengths(FS,VS));
dl=mean(patchEdgeLengths(FL,VL));
ligamentOrthogonalOffset=mean([ds dl])/2;
%ligamentOrthogonalOffset=mean([ds dl])/1;

VS_attach_curve=VS(indCurve_S,:)+ligamentOrthogonalOffset*NSv(indCurve_S,:);
VL_attach_curve=VL(indCurve_L,:)+ligamentOrthogonalOffset*NLv(indCurve_L,:);

VS_attach=[VS(indCurve_S,:); VS_attach_curve];
VL_attach=[VL(indCurve_L,:); VL_attach_curve];

FS_attach=[(1:1:numel(indCurve_S))' (1:1:numel(indCurve_S))'+numel(indCurve_S) [2:1:numel(indCurve_S) 1]'+numel(indCurve_S) [2:1:numel(indCurve_S) 1]'];
FS_attach=[FS_attach(:,[1 2 3]); FS_attach(:,[3 4 1]);];
FL_attach=[(1:1:numel(indCurve_L))' (1:1:numel(indCurve_L))'+numel(indCurve_L) [2:1:numel(indCurve_L) 1]'+numel(indCurve_L) [2:1:numel(indCurve_L) 1]'];
FL_attach=[FL_attach(:,[1 2 3]); FL_attach(:,[3 4 1]);];
FL_attach=fliplr(FL_attach); %Invert this face

d=sqrt(sum((mean(VL(indCurve_L,:),1)-mean(VS(indCurve_S,:),1)).^2));

cPar.closeLoopOpt=1;
cPar.patchType='tri_slash';
%VX=(VS(indCurve_S,:)+VL(indCurve_L,:))/2;
VX1=VS(indCurve_S,:)+nS*d/2;
VX2=VL(indCurve_L,:)+nL*d/2;
VX=(VX1+VX2)/2;

[Fn1,Vn1]=polyLoftLinear(VX,VS_attach_curve,cPar);
[Fn2,Vn2]=polyLoftLinear(VL_attach_curve,VX,cPar);

[F_ligament,V_ligament]=joinElementSets({Fn1,Fn2},{Vn1,Vn2});
[F_ligament,V_ligament]=mergeVertices(F_ligament,V_ligament);

Ebn=patchBoundary(F_ligament,V_ligament);
cParSmoothLoft.RigidConstraints=unique(Ebn(:));
[V_ligament]=patchSmooth(F_ligament,V_ligament,[],cParSmoothLoft);

%
[F_ligament,V_ligament]=joinElementSets({F_ligament,FS_attach,FL_attach},{V_ligament,VS_attach,VL_attach});
[F_ligament,V_ligament]=mergeVertices(F_ligament,V_ligament);

Ebn=patchBoundary(F_ligament,V_ligament);
cParSmoothTransition.RigidConstraints=unique(Ebn(:));
[V_ligament]=patchSmooth(F_ligament,V_ligament,[],cParSmoothTransition);

Ebn=patchBoundary(F_ligament,V_ligament);
cParSmoothTransition.RigidConstraints=unique(Ebn(:));
[V_ligament]=patchSmooth(F_ligament,V_ligament,[],cParSmoothTransition);


%%
numStepsCurve=50; %Number of steps for the curve
p1=mean(VL(indCurve_L,:),1); %First point
n1=nL; %First direction vector
p2=mean(VS(indCurve_S,:),1); %End point
n2=-nS; %End direction vector
f=0.05; %Extent of tangential nature to boundary curves, surface will remain approximately orthogonal to input curves for f*distance between curves
[Vg]=sweepCurveSmooth(p1,p2,n1,n2,numStepsCurve,csapsSmoothPar,f);
% [Fn,Vn,Cn]=sweepLoft(VL(indCurve_L,:),VS(indCurve_S,:),n1,n2,Vg,[],0,0);

%Visualization
cFigure; hold on;

gpatch(FL,VL,'w','none',0.5);
gpatch(FS,VS,'w','none',0.5);


plotV(VL(indCurve_L,:),'g.-','MarkerSize',25,'LineWidth',3);
plotV(VS(indCurve_S,:),'g.-','MarkerSize',25,'LineWidth',3);
plotV(Vg,'r.-','MarkerSize',25,'LineWidth',3);

gpatch(F_ligament,V_ligament,'gw','k',1);
patchNormPlot(F_ligament,V_ligament);
camlight headlight;
axisGeom(gca,fontSize);
view(viewSpec);
drawnow;

%% Fix overlap between bones and ligament
if fixOverlapOpt==1
    indBoundary=unique(patchBoundary(F_ligament,V_ligament));
    
    D=minDist(V_ligament,VL);
    logicClose=D<distanceTooClose;
    logicClose(indBoundary)=0;
    if any(logicClose)
        while 1
            [N,~,Nv]=patchNormal(F_ligament,V_ligament);
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
    
    
    %% Fix overlap
    
    [F_ligament,V_ligament]=fixOverlap(F_ligament,V_ligament,VL,distanceTooClose,stepSizeTooCloseFix);
    [F_ligament,V_ligament]=fixOverlap(F_ligament,V_ligament,VS,distanceTooClose,stepSizeTooCloseFix);
end

%%

cFigure; hold on;
gpatch(FL,VL,'w','none',0.5);
patchNormPlot(FL,VL);
gpatch(FS,VS,'w','none',0.5);
patchNormPlot(FS,VS);

plotV(VL(indCurve_L,:),'g.-','MarkerSize',25,'LineWidth',3);
plotV(VS(indCurve_S,:),'g.-','MarkerSize',25,'LineWidth',3);

plotV(Vg,'r.-','MarkerSize',25,'LineWidth',3);

gpatch(F_ligament,V_ligament,'gw','k',1);
patchNormPlot(F_ligament,V_ligament);
camlight headlight;
axisGeom(gca,fontSize);
view(viewSpec);
drawnow;

%Joint elements
[FT,VT,CT]=joinElementSets({F_ligament,FL(CL==1,:),FS(CS==1,:)},{V_ligament,VL,VS});
[FT,VT]=patchCleanUnused(FT,VT);
[FT,VT]=mergeVertices(FT,VT);
[V_regions]=getInnerPoint(FT,VT);

%%
cFigure; hold on;

gpatch(FL,VL,'w','none',0.5);
gpatch(FS,VS,'w','none',0.5);

gpatch(FT,VT,'g','k',1);
plotV(V_regions,'r.','MarkerSize',50);
camlight headlight;
axisGeom(gca,fontSize);
view(viewSpec);

%% Create solid mesh using tetgen

inputStruct.stringOpt='-pq1.2AaY';
inputStruct.Faces=FT;
inputStruct.Nodes=VT;
inputStruct.holePoints=[];
inputStruct.faceBoundaryMarker=CT; %Face boundary markers
inputStruct.regionPoints=V_regions; %region points
inputStruct.regionA=tetVolMeanEst(FT,VT); %Volume for regular tets
inputStruct.minRegionMarker=2; %Minimum region marker

% Mesh model using tetrahedral elements using tetGen 
[meshStruct]=runTetGen(inputStruct); %Run tetGen 

% Access model element and patch data
Fb=meshStruct.facesBoundary;
Cb=meshStruct.boundaryMarker;
V_tet=meshStruct.nodes;
CE=meshStruct.elementMaterialID;
E=meshStruct.elements;


%% Visualizing mesh using |meshView|, see also |anim8|

hFig=cFigure; 
hold on; 
title('Cut view of solid mesh','FontSize',fontSize);
gpatch(FL,VL,'kw','none',0.5);
gpatch(FS,VS,'kw','none',0.5);
gpatch(Fb,V_tet,'kw','none',0.5); 
optionStruct.hFig=hFig;
meshView(meshStruct,optionStruct);
axisGeom(gca,fontSize);
view(viewSpec);
drawnow;

%% Join and merge node sets
% Joining sets
V=[VL;VS;V_tet];

%Fixing indices for scaphoid
FS=FS+size(VL,1);

%Fixing indices for solid mesh and associated faces for plotting
Fb=Fb+size(VL,1)+size(VS,1);
E=E+size(VL,1)+size(VS,1);

% Merging
numDigitsMerge=6-numOrder(mean(patchEdgeLengths(Fb,V))); %base number of digits on mean
[~,indKeep,indFix]=unique(pround(V,numDigitsMerge),'rows');
V=V(indKeep,:);
FS=indFix(FS);
FL=indFix(FL);
Fb=indFix(Fb);
E=indFix(E);

F_tet=element2patch(E);

%Fibers
VE=patchCentre(E,V);
[~,indMin1]=minDist(VE,Vg);
indMin2=indMin1-1;
logicInvalid=indMin2<1;
indMin2(logicInvalid)=indMin1(logicInvalid)+1;
V_fib=Vg(indMin1,:)-Vg(indMin2,:);
V_fib(logicInvalid,:)=-V_fib(logicInvalid,:); 
V_fib=vecnormalize(V_fib);
[aFib,dFib]=vectorOrthogonalPair(V_fib);

%%
hFig=cFigure; 
hold on; 
%title('Cut view of solid mesh','FontSize',fontSize);
gpatch(FL,V,'w','none',0.5);
gpatch(FS,V,'w','none',0.5);
gpatch(Fb,V,'gw','none',0.5); 
plotV(Vg,'r.-','MarkerSize',25,'LineWidth',3);

[h]=quiverVec(VE,V_fib,0.1,'k');
[h]=quiverVec(VE,aFib,0.1,'r');
[h]=quiverVec(VE,dFib,0.1,'g');

% gpatch(F_tet,V,'kw','none',0.5); 
axisGeom(gca,fontSize);
view(viewSpec);
camlight headlight;
drawnow;

%%
hFig=cFigure; 
hold on; 
gpatch(Fb,V,'gw','none',0.5); 
plotV(Vg,'r.-','MarkerSize',25,'LineWidth',3);

[h]=quiverVec(VE,V_fib,0.1,'k');
%[h]=quiverVec(VE,aFib,0.1,'r');
%[h]=quiverVec(VE,dFib,0.1,'g');

axisGeom(gca,fontSize);
view(viewSpec);
camlight headlight;
drawnow;


%%
hFig=cFigure; 
hold on; 
gpatch(F_tet,V,'rw','none',0.5); 
axisGeom(gca,fontSize);
view(viewSpec);
camlight headlight;
drawnow;
%%



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