%% STUDY_04_wrist_scaffold_multiple_fiber_evaluation
% Below is a demonstration for:
% 
% * Finite element analysis of the performance of additively manufactured scaffolds for scapholunate ligament reconstruction
%___________________________________________________________________________________________________________________________
%  This study presents a patient-specific computational biomechanical evaluation of the effect of scaffold length, 
%  and positioning of the bone attachment sites. Through segmentation and image processing of medical image data 
%  for natural wrist motion, detailed 3D geometries as well as patient-specific physiological wrist motion could 
%  be derived. This data formed the input for detailed finite element analysis, enabling computational of scaffold 
%  stress and strain distributions, which are key predictors of scaffold structural integrity.
%___________________________________________________________________________________________________________________________
% * febio_spec version 3.0
% * febio, FEBio
% * hexahedral elements, hex8
% * static, solid
% * hyperelastic, Ogden
% * displacement logfile
% * stress logfile

clear; close all; clc;

% Plot settings
fontSize=40;
numSmoothStepsMain=1;
cParSmoothMain.n=numSmoothStepsMain;
cParSmoothMain.Method='HC';

%%
% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'2021_Nataliya_Perevoshchikova_wrist_project','data','temp');
loadPathSurface=fullfile(defaultFolder,'2021_Nataliya_Perevoshchikova_wrist_project','data','STL');
loadPathMotion=fullfile(defaultFolder,'2021_Nataliya_Perevoshchikova_wrist_project','data','motion');

%% Building a quadrilateral circular mesh
d=0.35;
h=7;
xSpacing=0.93;
ySpacing=0.525;
nCopies_x=4;
nCopies_y=7;
r=d/2;
ne=2; %Elements in radius
f=0.6; %Fraction (with respect to outer radius) where central square appears

pointSpacingHeight=0.1;
numStepsSweep=ceil(h/pointSpacingHeight);

xc=zeros(numStepsSweep,1);
yc=zeros(numStepsSweep,1);
zc=linspace(-h/2,h/2,numStepsSweep)';
Vc=[xc yc zc];

%SweepLoft
%Create the mesh
[Fq,Vq]=discQuadMesh(ne,r,f);

Vq(:,3)=0;
Vq1=Vq;
Vq1(:,3)=Vq1(:,3)-h/2;
Vq2=Vq;
Vq2(:,3)=Vq2(:,3)+h/2;

% Visualizing mesh
cFigure; hold on;
title('Single fiber','FontSize',fontSize);
plotV(Vq1,'k.','lineWidth',5,'MarkerSize',25);
plotV(Vq2,'r.','lineWidth',5,'MarkerSize',25);
plotV(Vc,'g.-','lineWidth',5,'MarkerSize',25);

axisGeom;
camlight headlight;
drawnow;

[~,~,~,S]=sweepLoft(Vq1,Vq2,[0 0 1],[0 0 1],Vc,numStepsSweep,0,0);

X=S.X'; Y=S.Y'; Z=S.Z'; %Coordinate matrices
V=[X(:) Y(:) Z(:)]; %Create node list

I=size(Vq,1)*((1:1:numStepsSweep)-1);
I=I(ones(size(Fq,1),1),:);
I=I(:);

FQ=repmat(Fq,numStepsSweep,1)+I(:,ones(size(Fq,2),1));
Efib=[FQ(1:end-size(Fq,1),:) FQ(size(Fq,1)+1:end,:)]; %The hexahedral elements

[Efib,V,ind1,ind2]=mergeVertices(Efib,V); %Merge nodes (start and end are not shared yet)
%%
nCopies=nCopies_x*nCopies_y;
VC=repmat({V},nCopies,1);
EC=repmat({Efib},nCopies,1);
numNodesBar=size(V,1);
[Efib,V]=joinElementSets(EC,VC);

c=1;

for qx=1:1:nCopies_x
    for qy=1:1:nCopies_y
        if c>1
            ind1=((c-1)*numNodesBar)+1;
            ind2=(c*numNodesBar);
            V(ind1:ind2,1)=V(ind1:ind2,1)+(qx-1)*xSpacing;
            V(ind1:ind2,2)=V(ind1:ind2,2)+(qy-1)*ySpacing;
        end
        c=c+1;
    end
end

[F]=element2patch(Efib);

%% Get boundary conditions faces
%Get boundary faces
ind=tesBoundary(F,V);
Fb=F(ind,:);
%Get tops/bottoms
Nb=patchNormal(Fb,V);

zVec=[0 0 1];
d=dot(Nb,zVec(ones(size(Nb,1),1),:),2);
Z=V(:,3);
ZF=mean(Z(Fb),2);
logicTop_Fb=(d>0.9) & ZF>=(max(V(:,3))-eps(1));
logicBottom_Fb=(d<-0.9) & ZF<=(min(V(:,3))+eps(1));
F_start=Fb(logicTop_Fb,:);
F_end=Fb(logicBottom_Fb,:);


logic_Lig=(ZF>-1.5) & ZF<=1.5;
F_Lig=Fb(logic_Lig,:);

%%
% Visualizing mesh
cFigure; hold on;
title('Multiple fibers','FontSize',fontSize);
patch('faces',Fb,'vertices',V,'FaceColor','g','FaceAlpha',0.3,'EdgeColor','None');
plotV(V(F_end',:),'k.','lineWidth',5,'MarkerSize',20);
plotV(V(F_start',:),'r.','lineWidth',5,'MarkerSize',20);

camlight headlight;
axisGeom(gca,fontSize);
drawnow;

F_L=F_end;
F_S=F_start;

%% Scaphoid and lunate at first frame from TriPlane, Brown University
loadName=fullfile(loadPathSurface,['modelGeometry.mat']);
dataStruct=load(loadName);

FS_triP=dataStruct.FS;
FL_triP=dataStruct.FL;
FC_triP=dataStruct.FC;
FR_triP=dataStruct.FR;

VL_triP=dataStruct.VL;
VS_triP=dataStruct.VS;
VC_triP=dataStruct.VC;
VR_triP=dataStruct.VR;

%%%
hFig=cFigure; 
hold on; 
title('Radiography','FontSize',fontSize);

gpatch(FS_triP,VS_triP,'w','none',0.5);
gpatch(FL_triP,VL_triP,'w','none',0.5);
gpatch(FC_triP,VC_triP,'w','none',0.5);
gpatch(FR_triP,VR_triP,'w','none',0.5);

grid on; axis equal;
axisGeom(gca,fontSize);
camlight headlight;
drawnow; 


%% Cutting Radius
pointSpacingR=mean(patchEdgeLengths(FR_triP,VR_triP));
%% Cut Capitate bone
% The capitate is cut in the x direction. 

%Create a logic for cutting away faces
max_X1=max(VR_triP(:,1))-70*pointSpacingR; %Max x-level used for cutting
logicVertices=VR_triP(:,1)<max_X1; %Logic for the points below this level
logicFaces=all(logicVertices(FR_triP),2); %Logic for the faces
logicFaces=triSurfLogicSharpFix(FR_triP,logicFaces,3); %Altered logic so it is smoother

%% 
% Visualize
cFigure; hold on;
gpatch(FS_triP,VS_triP,'w','none',0.5);
gpatch(FL_triP,VL_triP,'w','none',0.5);
gpatch(FR_triP,VR_triP,logicFaces,'none',0.5);
gpatch(FC_triP,VC_triP,'w','none',0.5);

camlight('headlight'); 
axisGeom(gca,fontSize);
drawnow; 

%%
% Cut away faces using logic
FR_triP=FR_triP(logicFaces,:); %The faces to keep
[FR_triP,VR_triP]=patchCleanUnused(FR_triP,VR_triP); %Remove unused points

%Attempt to self triangulate potentially jagged edge
Eb=patchBoundary(FR_triP,VR_triP); %Get boundary edges
indBoundary=edgeListToCurve(Eb); %Convert boundary edges to a curve list
indBoundary=indBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop
angleThreshold=pi*(120/180); %threshold for self triangulation 
[FR_triP,VR_triP,indBoundaryTop]=triSurfSelfTriangulateBoundary(FR_triP,VR_triP,indBoundary,angleThreshold,1);

%Force boundary to have the max X level chosen
VR_triP(indBoundaryTop,1)=max_X1;

%% 
% Visualize
cFigure; hold on;
gpatch(FS_triP,VS_triP,'w','none',0.5);
gpatch(FL_triP,VL_triP,'w','none',0.5);
gpatch(FC_triP,VC_triP,'w','none',0.5);
gpatch(FR_triP,VR_triP,'w','k',0.5);
plotV(VR_triP(indBoundaryTop,:),'r-','LineWidth',3);

camlight('headlight'); 
axisGeom(gca,fontSize);
drawnow; 

%% Close over top of radius
% The top boundary curve of the cut surface is filled with triangles. This
% is a 2D method. The x-coordinate is added after. 
[F2t,V2t]=regionTriMesh2D({VR_triP(indBoundaryTop,[2 3])},pointSpacingR,0,0); 
V2t(:,3)=V2t(:,2);
V2t(:,2)=V2t(:,1);
V2t(:,1)=max_X1; %Add/set X-level

%% 
% Visualize
cFigure; hold on;
gpatch(FS_triP,VS_triP,'w','none',0.5);
gpatch(FL_triP,VL_triP,'w','none',0.5);
gpatch(FC_triP,VC_triP,'w','none',0.5);
gpatch(FR_triP,VR_triP,'rw','k',0.5);
gpatch(F2t,V2t,'gw','k',1);
plotV(VR_triP(indBoundaryTop,:),'r-','LineWidth',3);
camlight('headlight'); 
axisGeom(gca,fontSize);
drawnow; 

%% Joining surface features
% Add all surface sets together in joint list of faces, vertices
[FR_triP,VR_triP]=joinElementSets({FR_triP,F2t},{VR_triP,V2t});
%% Merge shared nodes
% The join operation only adds the sets together. Nodes with the same
% coordinates are not seen as the same yet and need to be merged. 
[FR_triP,VR_triP]=mergeVertices(FR_triP,VR_triP); %Merge nodes

%Remove tri connected
[FS,VS]=triSurfRemoveThreeConnect(FS_triP,VS_triP);
[FL,VL]=triSurfRemoveThreeConnect(FL_triP,VL_triP);
[FC,VC]=triSurfRemoveThreeConnect(FC_triP,VC_triP);
[FR,VR]=triSurfRemoveThreeConnect(FR_triP,VR_triP);
%Smoothen
[VL]=patchSmooth(FL,VL,[],cParSmoothMain);
[VS]=patchSmooth(FS,VS,[],cParSmoothMain);
[VR]=patchSmooth(FR,VR,[],cParSmoothMain);
[VC]=patchSmooth(FC,VC,[],cParSmoothMain);
%Refine
[FS,VS]=subtri(FS,VS,1);
[FL,VL]=subtri(FL,VL,1);
[FR,VR]=subtri(FR,VR,1);
[FC,VC]=subtri(FC,VC,1);
%Smoothen
[VL]=patchSmooth(FL,VL,[],cParSmoothMain);
[VS]=patchSmooth(FS,VS,[],cParSmoothMain);
[VR]=patchSmooth(FR,VR,[],cParSmoothMain);
[VC]=patchSmooth(FC,VC,[],cParSmoothMain);
%%%%

%Construct placement according to natural SLIL positioning
alphaRot=pi/2;
R=euler2DCM([alphaRot 0 0]); %The rotation tensor for each step
V=V*R; %Rotated 1 part for visualization of stepwise amount
 
alphaRot=pi/2;
R=euler2DCM([0 alphaRot 0]); %The rotation tensor for each step
V=V*R; %Rotated 1 part for visualization of stepwise amount

alphaRot=-pi/15;
R=euler2DCM([0 0 alphaRot]); %The rotation tensor for each step
V=V*R; %Rotated 1 part for visualization of stepwise amount

V(:,3)=V(:,3)-1.5;
V(:,1)=V(:,1)-15;
V(:,2)=V(:,2)-3.2;

hFig=cFigure; 
hold on; 

gpatch(FS,VS,'w','none',0.7);
gpatch(FL,VL,'w','none',0.7);

gpatch(FC,VC,'w','none',0.5);
gpatch(FR,VR,'w','none',0.5);

gpatch(Fb,V,'gw','k',1);

grid on; axis equal;
axisGeom(gca,fontSize);
camlight headlight;
drawnow;    

%% Joining node sets
[E_F,V,C]=joinElementSets({Efib,FS,FL,FC,FR},{V,VS,VL,VC,VR});
E=E_F{1};
FS1=E_F{2};
FL1=E_F{3};
FC1=E_F{4};
FR1=E_F{5};

% hFig=cFigure; 
% hold on; 
% 
% gpatch(FS1,V,'w','none',0.7);
% gpatch(FL1,V,'w','none',0.7);
% 
% gpatch(FC1,V,'w','none',0.5);
% gpatch(FR1,V,'w','none',0.5);
% 
% gpatch(E,V,'gw','k',1);
% 
% grid on; axis equal;
% axisGeom(gca,fontSize);
% camlight headlight;
% drawnow;  


%% FEA control settings
Step=5;
numTimeSteps=round((1989-1000)/Step);%Number of time steps desired
max_refs=300; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=150; %Optimum number of iterations
max_retries=50; %Maximum number of retires
dtmin=(1/numTimeSteps)/10; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
tStep=0.005*Step; %Duration of step
min_residual=1e-20;
Mag_val=1;


%% Control parameters
% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='3.0'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

%Control section
febio_spec.Control.analysis='STATIC';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=tStep;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.solver.max_ups=max_ups;
febio_spec.Control.solver.min_residual=min_residual;


%Material section
materialName1='Material1';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.ATTR.name=materialName1;

%Material section
C1=57.21894783089972;
C2=57.23071414361289;
K=11444.966197451262;
M1=2.00012552045284;
M2=-1.9997730268725347;
t1=25.365;
g1=1.23339;
        
materialOption='elastic';%'elastic_viscoelastic'
switch materialOption
    case 'elastic'
        febio_spec.Material.material{1}.ATTR.type='Ogden';        
        febio_spec.Material.material{1}.c1=C1;
        febio_spec.Material.material{1}.m1=M1;
        febio_spec.Material.material{1}.c2=C2;
        febio_spec.Material.material{1}.m2=M2;
        febio_spec.Material.material{1}.k=K;
    case 'elastic_viscoelastic'
        %Viscoelastic part
        febio_spec.Material.material{1}.ATTR.type='uncoupled viscoelastic';
        febio_spec.Material.material{1}.g1=g1;
        febio_spec.Material.material{1}.t1=t1;

        %Elastic part
        febio_spec.Material.material{1}.elastic{1}.ATTR.type='Ogden';
        febio_spec.Material.material{1}.elastic{1}.c1=C1;
        febio_spec.Material.material{1}.elastic{1}.m1=M1;
        febio_spec.Material.material{1}.elastic{1}.c2=C2;
        febio_spec.Material.material{1}.elastic{1}.m2=M2;
        febio_spec.Material.material{1}.elastic{1}.k=K;
        
        
end

materialName2='Lunate_side';
febio_spec.Material.material{2}.ATTR.name=materialName2;
febio_spec.Material.material{2}.ATTR.type='rigid body';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.density=1;
febio_spec.Material.material{2}.center_of_mass=[0 0 0];

materialName3='Scaphoid_side';
febio_spec.Material.material{3}.ATTR.name=materialName3;
febio_spec.Material.material{3}.ATTR.type='rigid body';
febio_spec.Material.material{3}.ATTR.id=3;
febio_spec.Material.material{3}.density=1;
febio_spec.Material.material{3}.center_of_mass=[0 0 0];

materialName4='Lunate';
febio_spec.Material.material{4}.ATTR.name=materialName4;
febio_spec.Material.material{4}.ATTR.type='rigid body';
febio_spec.Material.material{4}.ATTR.id=4;
febio_spec.Material.material{4}.density=1;
febio_spec.Material.material{4}.center_of_mass=[0 0 0];

materialName5='Scaphoid';
febio_spec.Material.material{5}.ATTR.name=materialName5;
febio_spec.Material.material{5}.ATTR.type='rigid body';
febio_spec.Material.material{5}.ATTR.id=5;
febio_spec.Material.material{5}.density=1;
febio_spec.Material.material{5}.center_of_mass=[0 0 0];

materialName6='Capitate';
febio_spec.Material.material{6}.ATTR.name=materialName6;
febio_spec.Material.material{6}.ATTR.type='rigid body';
febio_spec.Material.material{6}.ATTR.id=6;
febio_spec.Material.material{6}.density=1;
febio_spec.Material.material{6}.center_of_mass=[0 0 0];

materialName7='Radius';
febio_spec.Material.material{7}.ATTR.name=materialName7;
febio_spec.Material.material{7}.ATTR.type='rigid body';
febio_spec.Material.material{7}.ATTR.id=7;
febio_spec.Material.material{7}.density=1;
febio_spec.Material.material{7}.center_of_mass=[0 0 0];

%Geometry section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% % -> Elements
partName1='Scaffold';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of the element set
febio_spec.Mesh.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type of this set
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E;

partName2='Lunate side';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of the element set
febio_spec.Mesh.Elements{2}.ATTR.mat=2; %material index for this set
febio_spec.Mesh.Elements{2}.ATTR.type='quad4'; %Element type of this set
febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E,1)+(1:1:size(F_L,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=F_L;

partName3='Scaphoid side';
febio_spec.Mesh.Elements{3}.ATTR.name=partName3; %Name of the element set
febio_spec.Mesh.Elements{3}.ATTR.mat=3; %material index for this set
febio_spec.Mesh.Elements{3}.ATTR.type='quad4'; %Element type of this set
febio_spec.Mesh.Elements{3}.elem.ATTR.id=size(E,1)+(1:1:size(F_S,1))'; %Element id's
febio_spec.Mesh.Elements{3}.elem.VAL=F_S;

partName4='Lunate';
febio_spec.Mesh.Elements{4}.ATTR.name=partName4; %Name of the element set
febio_spec.Mesh.Elements{4}.ATTR.mat=4; %material index for this set
febio_spec.Mesh.Elements{4}.ATTR.type='tri3'; %Element type of this set
febio_spec.Mesh.Elements{4}.elem.ATTR.id=size(E,1)+(1:1:size(FL1,1))'; %Element id's
febio_spec.Mesh.Elements{4}.elem.VAL=FL1;


partName5='Scaphoid';
febio_spec.Mesh.Elements{5}.ATTR.name=partName5; %Name of the element set
febio_spec.Mesh.Elements{5}.ATTR.mat=5; %material index for this set
febio_spec.Mesh.Elements{5}.ATTR.type='tri3'; %Element type of this set
febio_spec.Mesh.Elements{5}.elem.ATTR.id=size(E,1)+(1:1:size(FS1,1))'; %Element id's
febio_spec.Mesh.Elements{5}.elem.VAL=FS1;

%Capitate
partName6='Capitate';
febio_spec.Mesh.Elements{6}.ATTR.name=partName6; %Name of the element set
febio_spec.Mesh.Elements{6}.ATTR.mat=6; %material index for this set
febio_spec.Mesh.Elements{6}.ATTR.type='tri3'; %Element type of this set
febio_spec.Mesh.Elements{6}.elem.ATTR.id=size(E,1)+(1:1:size(FC1,1))'; %Element id's
febio_spec.Mesh.Elements{6}.elem.VAL=FC1;

%Radius
partName7='Radius';
febio_spec.Mesh.Elements{7}.ATTR.name=partName7; %Name of the element set
febio_spec.Mesh.Elements{7}.ATTR.mat=7; %material index for this set
febio_spec.Mesh.Elements{7}.ATTR.type='tri3'; %Element type of this set
febio_spec.Mesh.Elements{7}.elem.ATTR.id=size(E,1)+(1:1:size(FR1,1))'; %Element id's
febio_spec.Mesh.Elements{7}.elem.VAL=FR1;

%MeshDomains section
febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

febio_spec.MeshDomains.ShellDomain{1}.ATTR.name=partName2;
febio_spec.MeshDomains.ShellDomain{1}.ATTR.mat=materialName2;

febio_spec.MeshDomains.ShellDomain{2}.ATTR.name=partName3;
febio_spec.MeshDomains.ShellDomain{2}.ATTR.mat=materialName3;

febio_spec.MeshDomains.ShellDomain{3}.ATTR.name=partName4;
febio_spec.MeshDomains.ShellDomain{3}.ATTR.mat=materialName4;

febio_spec.MeshDomains.ShellDomain{4}.ATTR.name=partName5;
febio_spec.MeshDomains.ShellDomain{4}.ATTR.mat=materialName5;

febio_spec.MeshDomains.ShellDomain{5}.ATTR.name=partName6;
febio_spec.MeshDomains.ShellDomain{5}.ATTR.mat=materialName6;

febio_spec.MeshDomains.ShellDomain{6}.ATTR.name=partName7;
febio_spec.MeshDomains.ShellDomain{6}.ATTR.mat=materialName7;

% Boundary conditions
i=1;
for count=1:numTimeSteps

    FileData = load(fullfile(loadPathMotion, sprintf('Flex_Ext_%d_TriPlane_Lunate.mat',count)));
    MatT_Lunate = FileData.Trans_L;
    
    FileData = load(fullfile(loadPathMotion, sprintf('Flex_Ext_%d_TriPlane_Scaphoid.mat',count)));
    MatT_Scaphoid = FileData.Trans_S;
    
    FileData = load(fullfile(loadPathMotion, sprintf('Flex_Ext_%d_TriPlane_Capitate.mat',count)));
    MatT_Capitate = FileData.Trans_C;

    
    RLunate                       = MatT_Lunate(1:3,1:3);%radians
    RScaphoid                     = MatT_Scaphoid(1:3,1:3); %radians
    RCapitate                     = MatT_Capitate(1:3,1:3); %radians
    
    Prescribe_Lunate{i}       = MatT_Lunate(1:3,end);
    Prescribe_Scaphoid{i}     = MatT_Scaphoid(1:3,end);
    Prescribe_Capitate{i}     = MatT_Capitate(1:3,end);

    axangL                    = rotm2axang(RLunate);%Axis-Angle {[x, y, z], angle (radians)}
    AngleLunate{i}            = axangL(1:3)*axangL(4);%angle *rot_axis; Axis with angle magnitude (radians) [x, y, z]
    axangS                    = rotm2axang(RScaphoid);
    AngleScaphoid{i}          = axangS(1:3)*axangS(4);%angle*rot_axis;
    axangC                    = rotm2axang(RCapitate);
    AngleCapitate{i}          = axangC(1:3)*axangC(4);%angle*rot_axis;    
    tStartNow(i)=count*tStep;
    i=i+1;
end

%Boundary conditions
%Rigid section 
% %Start and End
febio_spec.Rigid.rigid_constraint{1}.ATTR.name='Rigid_Lunate_X';
febio_spec.Rigid.rigid_constraint{1}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{1}.rb=2;
febio_spec.Rigid.rigid_constraint{1}.dof='Rx';
febio_spec.Rigid.rigid_constraint{1}.value.ATTR.lc=1;
febio_spec.Rigid.rigid_constraint{1}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{1}.relative=0;

febio_spec.Rigid.rigid_constraint{2}.ATTR.name='Rigid_Lunate_Y';
febio_spec.Rigid.rigid_constraint{2}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{2}.rb=2;
febio_spec.Rigid.rigid_constraint{2}.dof='Ry';
febio_spec.Rigid.rigid_constraint{2}.value.ATTR.lc=2;
febio_spec.Rigid.rigid_constraint{2}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{2}.relative=0;

febio_spec.Rigid.rigid_constraint{3}.ATTR.name='Rigid_Lunate_Z';
febio_spec.Rigid.rigid_constraint{3}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{3}.rb=2;
febio_spec.Rigid.rigid_constraint{3}.dof='Rz';
febio_spec.Rigid.rigid_constraint{3}.value.ATTR.lc=3;
febio_spec.Rigid.rigid_constraint{3}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{3}.relative=0;

febio_spec.Rigid.rigid_constraint{4}.ATTR.name='Rigid_Lunate_Ru';
febio_spec.Rigid.rigid_constraint{4}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{4}.rb=2;
febio_spec.Rigid.rigid_constraint{4}.dof='Ru';
febio_spec.Rigid.rigid_constraint{4}.value.ATTR.lc=4;
febio_spec.Rigid.rigid_constraint{4}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{4}.relative=0;

febio_spec.Rigid.rigid_constraint{5}.ATTR.name='Rigid_Lunate_Rv';
febio_spec.Rigid.rigid_constraint{5}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{5}.rb=2;
febio_spec.Rigid.rigid_constraint{5}.dof='Rv';
febio_spec.Rigid.rigid_constraint{5}.value.ATTR.lc=5;
febio_spec.Rigid.rigid_constraint{5}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{5}.relative=0;

febio_spec.Rigid.rigid_constraint{6}.ATTR.name='Rigid_Lunate_Rw';
febio_spec.Rigid.rigid_constraint{6}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{6}.rb=2;
febio_spec.Rigid.rigid_constraint{6}.dof='Rw';
febio_spec.Rigid.rigid_constraint{6}.value.ATTR.lc=6;
febio_spec.Rigid.rigid_constraint{6}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{6}.relative=0;

febio_spec.Rigid.rigid_constraint{7}.ATTR.name='Rigid_Scaphoid_X';
febio_spec.Rigid.rigid_constraint{7}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{7}.rb=3;
febio_spec.Rigid.rigid_constraint{7}.dof='Rx';
febio_spec.Rigid.rigid_constraint{7}.value.ATTR.lc=7;
febio_spec.Rigid.rigid_constraint{7}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{7}.relative=0;

febio_spec.Rigid.rigid_constraint{8}.ATTR.name='Rigid_Scaphoid_Y';
febio_spec.Rigid.rigid_constraint{8}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{8}.rb=3;
febio_spec.Rigid.rigid_constraint{8}.dof='Ry';
febio_spec.Rigid.rigid_constraint{8}.value.ATTR.lc=8;
febio_spec.Rigid.rigid_constraint{8}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{8}.relative=0;

febio_spec.Rigid.rigid_constraint{9}.ATTR.name='Rigid_Scaphoid_Z';
febio_spec.Rigid.rigid_constraint{9}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{9}.rb=3;
febio_spec.Rigid.rigid_constraint{9}.dof='Rz';
febio_spec.Rigid.rigid_constraint{9}.value.ATTR.lc=9;
febio_spec.Rigid.rigid_constraint{9}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{9}.relative=0;

febio_spec.Rigid.rigid_constraint{10}.ATTR.name='Rigid_Scaphoid_RX';
febio_spec.Rigid.rigid_constraint{10}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{10}.rb=3;
febio_spec.Rigid.rigid_constraint{10}.dof='Ru';
febio_spec.Rigid.rigid_constraint{10}.value.ATTR.lc=10;
febio_spec.Rigid.rigid_constraint{10}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{10}.relative=0;

febio_spec.Rigid.rigid_constraint{11}.ATTR.name='Rigid_Scaphoid_RY';
febio_spec.Rigid.rigid_constraint{11}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{11}.rb=3;
febio_spec.Rigid.rigid_constraint{11}.dof='Rv';
febio_spec.Rigid.rigid_constraint{11}.value.ATTR.lc=11;
febio_spec.Rigid.rigid_constraint{11}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{11}.relative=0;

febio_spec.Rigid.rigid_constraint{12}.ATTR.name='Rigid_Scaphoid_RZ';
febio_spec.Rigid.rigid_constraint{12}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{12}.rb=3;
febio_spec.Rigid.rigid_constraint{12}.dof='Rw';
febio_spec.Rigid.rigid_constraint{12}.value.ATTR.lc=12;
febio_spec.Rigid.rigid_constraint{12}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{12}.relative=0;

% % -> Prescribed boundary conditions on the rigid body bones
febio_spec.Rigid.rigid_constraint{13}.ATTR.name='Rigid_Lunate_side_X';
febio_spec.Rigid.rigid_constraint{13}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{13}.rb=4;
febio_spec.Rigid.rigid_constraint{13}.dof='Rx';
febio_spec.Rigid.rigid_constraint{13}.value.ATTR.lc=1;
febio_spec.Rigid.rigid_constraint{13}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{13}.relative=0;

febio_spec.Rigid.rigid_constraint{14}.ATTR.name='Rigid_Lunate_side_Y';
febio_spec.Rigid.rigid_constraint{14}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{14}.rb=4;
febio_spec.Rigid.rigid_constraint{14}.dof='Ry';
febio_spec.Rigid.rigid_constraint{14}.value.ATTR.lc=2;
febio_spec.Rigid.rigid_constraint{14}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{14}.relative=0;

febio_spec.Rigid.rigid_constraint{15}.ATTR.name='Rigid_Lunate_side_Z';
febio_spec.Rigid.rigid_constraint{15}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{15}.rb=4;
febio_spec.Rigid.rigid_constraint{15}.dof='Rz';
febio_spec.Rigid.rigid_constraint{15}.value.ATTR.lc=3;
febio_spec.Rigid.rigid_constraint{15}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{15}.relative=0;

febio_spec.Rigid.rigid_constraint{16}.ATTR.name='Rigid_Lunate_side_RX';
febio_spec.Rigid.rigid_constraint{16}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{16}.rb=4;
febio_spec.Rigid.rigid_constraint{16}.dof='Ru';
febio_spec.Rigid.rigid_constraint{16}.value.ATTR.lc=4;
febio_spec.Rigid.rigid_constraint{16}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{16}.relative=0;

febio_spec.Rigid.rigid_constraint{17}.ATTR.name='Rigid_Lunate_side_RY';
febio_spec.Rigid.rigid_constraint{17}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{17}.rb=4;
febio_spec.Rigid.rigid_constraint{17}.dof='Rv';
febio_spec.Rigid.rigid_constraint{17}.value.ATTR.lc=5;
febio_spec.Rigid.rigid_constraint{17}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{17}.relative=0;

febio_spec.Rigid.rigid_constraint{18}.ATTR.name='Rigid_Lunate_side_RZ';
febio_spec.Rigid.rigid_constraint{18}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{18}.rb=4;
febio_spec.Rigid.rigid_constraint{18}.dof='Rw';
febio_spec.Rigid.rigid_constraint{18}.value.ATTR.lc=6;
febio_spec.Rigid.rigid_constraint{18}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{18}.relative=0;

% %Scaphoid
febio_spec.Rigid.rigid_constraint{19}.ATTR.name='Rigid_Scaphoid_side_X';
febio_spec.Rigid.rigid_constraint{19}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{19}.rb=5;
febio_spec.Rigid.rigid_constraint{19}.dof='Rx';
febio_spec.Rigid.rigid_constraint{19}.value.ATTR.lc=7;
febio_spec.Rigid.rigid_constraint{19}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{19}.relative=0;

febio_spec.Rigid.rigid_constraint{20}.ATTR.name='Rigid_Scaphoid_side_Y';
febio_spec.Rigid.rigid_constraint{20}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{20}.rb=5;
febio_spec.Rigid.rigid_constraint{20}.dof='Ry';
febio_spec.Rigid.rigid_constraint{20}.value.ATTR.lc=8;
febio_spec.Rigid.rigid_constraint{20}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{20}.relative=0;

febio_spec.Rigid.rigid_constraint{21}.ATTR.name='Rigid_Scaphoid_side_Z';
febio_spec.Rigid.rigid_constraint{21}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{21}.rb=5;
febio_spec.Rigid.rigid_constraint{21}.dof='Rz';
febio_spec.Rigid.rigid_constraint{21}.value.ATTR.lc=9;
febio_spec.Rigid.rigid_constraint{21}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{21}.relative=0;

febio_spec.Rigid.rigid_constraint{22}.ATTR.name='Rigid_Scaphoid_side_RX';
febio_spec.Rigid.rigid_constraint{22}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{22}.rb=5;
febio_spec.Rigid.rigid_constraint{22}.dof='Ru';
febio_spec.Rigid.rigid_constraint{22}.value.ATTR.lc=10;
febio_spec.Rigid.rigid_constraint{22}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{22}.relative=0;

febio_spec.Rigid.rigid_constraint{23}.ATTR.name='Rigid_Scaphoid_side_RY';
febio_spec.Rigid.rigid_constraint{23}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{23}.rb=5;
febio_spec.Rigid.rigid_constraint{23}.dof='Rv';
febio_spec.Rigid.rigid_constraint{23}.value.ATTR.lc=11;
febio_spec.Rigid.rigid_constraint{23}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{23}.relative=0;

febio_spec.Rigid.rigid_constraint{24}.ATTR.name='Rigid_Scaphoid_side_RZ';
febio_spec.Rigid.rigid_constraint{24}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{24}.rb=5;
febio_spec.Rigid.rigid_constraint{24}.dof='Rw';
febio_spec.Rigid.rigid_constraint{24}.value.ATTR.lc=12;
febio_spec.Rigid.rigid_constraint{24}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{24}.relative=0;

%%Capitate
febio_spec.Rigid.rigid_constraint{25}.ATTR.name='Rigid_Capitate_X';
febio_spec.Rigid.rigid_constraint{25}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{25}.rb=6;
febio_spec.Rigid.rigid_constraint{25}.dof='Rx';
febio_spec.Rigid.rigid_constraint{25}.value.ATTR.lc=13;
febio_spec.Rigid.rigid_constraint{25}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{25}.relative=0;

febio_spec.Rigid.rigid_constraint{26}.ATTR.name='Rigid_Capitate_Y';
febio_spec.Rigid.rigid_constraint{26}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{26}.rb=6;
febio_spec.Rigid.rigid_constraint{26}.dof='Ry';
febio_spec.Rigid.rigid_constraint{26}.value.ATTR.lc=14;
febio_spec.Rigid.rigid_constraint{26}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{26}.relative=0;

febio_spec.Rigid.rigid_constraint{27}.ATTR.name='Rigid_Capitate_Z';
febio_spec.Rigid.rigid_constraint{27}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{27}.rb=6;
febio_spec.Rigid.rigid_constraint{27}.dof='Rz';
febio_spec.Rigid.rigid_constraint{27}.value.ATTR.lc=15;
febio_spec.Rigid.rigid_constraint{27}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{27}.relative=0;

febio_spec.Rigid.rigid_constraint{28}.ATTR.name='Rigid_Capitate_RX';
febio_spec.Rigid.rigid_constraint{28}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{28}.rb=6;
febio_spec.Rigid.rigid_constraint{28}.dof='Ru';
febio_spec.Rigid.rigid_constraint{28}.value.ATTR.lc=16;
febio_spec.Rigid.rigid_constraint{28}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{28}.relative=0;

febio_spec.Rigid.rigid_constraint{29}.ATTR.name='Rigid_Capitate_RY';
febio_spec.Rigid.rigid_constraint{29}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{29}.rb=6;
febio_spec.Rigid.rigid_constraint{29}.dof='Rv';
febio_spec.Rigid.rigid_constraint{29}.value.ATTR.lc=17;
febio_spec.Rigid.rigid_constraint{29}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{29}.relative=0;

febio_spec.Rigid.rigid_constraint{30}.ATTR.name='Rigid_Capitate_RZ';
febio_spec.Rigid.rigid_constraint{30}.ATTR.type='prescribe';
febio_spec.Rigid.rigid_constraint{30}.rb=6;
febio_spec.Rigid.rigid_constraint{30}.dof='Rw';
febio_spec.Rigid.rigid_constraint{30}.value.ATTR.lc=18;
febio_spec.Rigid.rigid_constraint{30}.value.VAL=Mag_val;
febio_spec.Rigid.rigid_constraint{30}.relative=0;

% %Radius
febio_spec.Rigid.rigid_constraint{31}.ATTR.name='Rigid_Radiuas';
febio_spec.Rigid.rigid_constraint{31}.ATTR.type='fix';
febio_spec.Rigid.rigid_constraint{31}.rb=7;
febio_spec.Rigid.rigid_constraint{31}.dofs='Rx,Ry,Rz,Ru,Rv,Rw';

k=1;
for i=1:numTimeSteps+1
    
    if k==1
        Mat_X_Lunate(k,1)=0;
        Mat_X_Lunate(k,2)=0;
      
        Mat_Y_Lunate(k,1)=0;
        Mat_Y_Lunate(k,2)=0;
        
        Mat_Z_Lunate(k,1)=0;
        Mat_Z_Lunate(k,2)=0;
        
        Mat_RX_Lunate(k,1)=0;
        Mat_RX_Lunate(k,2)=0;
        
        Mat_RY_Lunate(k,1)=0;
        Mat_RY_Lunate(k,2)=0;
        
        Mat_RZ_Lunate(k,1)=0;
        Mat_RZ_Lunate(k,2)=0;
        
        Mat_X_Scaphoid(k,1)=0;
        Mat_X_Scaphoid(k,2)=0;
        
        Mat_Y_Scaphoid(k,1)=0;
        Mat_Y_Scaphoid(k,2)=0;
        
        Mat_Z_Scaphoid(k,1)=0;
        Mat_Z_Scaphoid(k,2)=0;
        
        Mat_RX_Scaphoid(k,1)=0;
        Mat_RX_Scaphoid(k,2)=0;
        
        Mat_RY_Scaphoid(k,1)=0;
        Mat_RY_Scaphoid(k,2)=0;
        
        Mat_RZ_Scaphoid(k,1)=0;
        Mat_RZ_Scaphoid(k,2)=0;
        
        Mat_X_Capitate(k,1)=0;
        Mat_X_Capitate(k,2)=0;
        
        Mat_Y_Capitate(k,1)=0;
        Mat_Y_Capitate(k,2)=0;
        
        Mat_Z_Capitate(k,1)=0;
        Mat_Z_Capitate(k,2)=0;
        
        Mat_RX_Capitate(k,1)=0;
        Mat_RX_Capitate(k,2)=0;
        
        Mat_RY_Capitate(k,1)=0;
        Mat_RY_Capitate(k,2)=0;
        
        Mat_RZ_Capitate(k,1)=0;
        Mat_RZ_Capitate(k,2)=0;
    else
        Mat_X_Lunate(k,1)=tStartNow(k-1);
        Mat_X_Lunate(k,2)=Prescribe_Lunate{k-1}(1);
        
        Mat_Y_Lunate(k,1)=tStartNow(k-1);
        Mat_Y_Lunate(k,2)=Prescribe_Lunate{k-1}(2);
        
        Mat_Z_Lunate(k,1)=tStartNow(k-1);
        Mat_Z_Lunate(k,2)=Prescribe_Lunate{k-1}(3);
        
        Mat_RX_Lunate(k,1)=tStartNow(k-1);
        Mat_RX_Lunate(k,2)=AngleLunate{k-1}(1);
        
        Mat_RY_Lunate(k,1)=tStartNow(k-1);
        Mat_RY_Lunate(k,2)=AngleLunate{k-1}(2);
        
        Mat_RZ_Lunate(k,1)=tStartNow(k-1);
        Mat_RZ_Lunate(k,2)=AngleLunate{k-1}(3);
        
        Mat_X_Scaphoid(k,1)=tStartNow(k-1);
        Mat_X_Scaphoid(k,2)=Prescribe_Scaphoid{k-1}(1);
        
        Mat_Y_Scaphoid(k,1)=tStartNow(k-1);
        Mat_Y_Scaphoid(k,2)=Prescribe_Scaphoid{k-1}(2);
        
        Mat_Z_Scaphoid(k,1)=tStartNow(k-1);
        Mat_Z_Scaphoid(k,2)=Prescribe_Scaphoid{k-1}(3);
        
        Mat_RX_Scaphoid(k,1)=tStartNow(k-1);
        Mat_RX_Scaphoid(k,2)=AngleScaphoid{k-1}(1);
        
        Mat_RY_Scaphoid(k,1)=tStartNow(k-1);
        Mat_RY_Scaphoid(k,2)=AngleScaphoid{k-1}(2);
        
        Mat_RZ_Scaphoid(k,1)=tStartNow(k-1);
        Mat_RZ_Scaphoid(k,2)=AngleScaphoid{k-1}(3);
        
        Mat_X_Capitate(k,1)=tStartNow(k-1);
        Mat_X_Capitate(k,2)=Prescribe_Capitate{k-1}(1);
        
        Mat_Y_Capitate(k,1)=tStartNow(k-1);
        Mat_Y_Capitate(k,2)=Prescribe_Capitate{k-1}(2);
        
        Mat_Z_Capitate(k,1)=tStartNow(k-1);
        Mat_Z_Capitate(k,2)=Prescribe_Capitate{k-1}(3);
        
        Mat_RX_Capitate(k,1)=tStartNow(k-1);
        Mat_RX_Capitate(k,2)=AngleCapitate{k-1}(1);
        
        Mat_RY_Capitate(k,1)=tStartNow(k-1);
        Mat_RY_Capitate(k,2)=AngleCapitate{k-1}(2);
        
        Mat_RZ_Capitate(k,1)=tStartNow(k-1);
        Mat_RZ_Capitate(k,2)=AngleCapitate{k-1}(3);
        
    end
   k=k+1; 
end

%Lunate
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{1}.points.point.VAL=Mat_X_Lunate;

febio_spec.LoadData.load_controller{2}.ATTR.id=2;
febio_spec.LoadData.load_controller{2}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{2}.points.point.VAL=Mat_Y_Lunate;

febio_spec.LoadData.load_controller{3}.ATTR.id=3;
febio_spec.LoadData.load_controller{3}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{3}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{3}.points.point.VAL=Mat_Z_Lunate;

febio_spec.LoadData.load_controller{4}.ATTR.id=4;
febio_spec.LoadData.load_controller{4}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{4}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{4}.points.point.VAL=Mat_RX_Lunate;

febio_spec.LoadData.load_controller{5}.ATTR.id=5;
febio_spec.LoadData.load_controller{5}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{5}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{5}.points.point.VAL=Mat_RY_Lunate;

febio_spec.LoadData.load_controller{6}.ATTR.id=6;
febio_spec.LoadData.load_controller{6}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{6}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{6}.points.point.VAL=Mat_RZ_Lunate;

% %Scaphoid   
febio_spec.LoadData.load_controller{7}.ATTR.id=7;
febio_spec.LoadData.load_controller{7}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{7}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{7}.points.point.VAL=Mat_X_Scaphoid;

febio_spec.LoadData.load_controller{8}.ATTR.id=8;
febio_spec.LoadData.load_controller{8}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{8}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{8}.points.point.VAL=Mat_Y_Scaphoid;

febio_spec.LoadData.load_controller{9}.ATTR.id=9;
febio_spec.LoadData.load_controller{9}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{9}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{9}.points.point.VAL=Mat_Z_Scaphoid;

febio_spec.LoadData.load_controller{10}.ATTR.id=10;
febio_spec.LoadData.load_controller{10}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{10}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{10}.points.point.VAL=Mat_RX_Scaphoid;

febio_spec.LoadData.load_controller{11}.ATTR.id=11;
febio_spec.LoadData.load_controller{11}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{11}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{11}.points.point.VAL=Mat_RY_Scaphoid;
       
febio_spec.LoadData.load_controller{12}.ATTR.id=12;
febio_spec.LoadData.load_controller{12}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{12}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{12}.points.point.VAL=Mat_RZ_Scaphoid;

% %Capitate
febio_spec.LoadData.load_controller{13}.ATTR.id=13;
febio_spec.LoadData.load_controller{13}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{13}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{13}.points.point.VAL=Mat_X_Capitate;

febio_spec.LoadData.load_controller{14}.ATTR.id=14;
febio_spec.LoadData.load_controller{14}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{14}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{14}.points.point.VAL=Mat_Y_Capitate;

febio_spec.LoadData.load_controller{15}.ATTR.id=15;
febio_spec.LoadData.load_controller{15}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{15}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{15}.points.point.VAL=Mat_Z_Capitate;

febio_spec.LoadData.load_controller{16}.ATTR.id=16;
febio_spec.LoadData.load_controller{16}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{16}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{16}.points.point.VAL=Mat_RX_Capitate;

febio_spec.LoadData.load_controller{17}.ATTR.id=17;
febio_spec.LoadData.load_controller{17}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{17}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{17}.points.point.VAL=Mat_RY_Capitate;
 
febio_spec.LoadData.load_controller{18}.ATTR.id=18;
febio_spec.LoadData.load_controller{18}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{18}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{18}.points.point.VAL=Mat_RZ_Capitate;

febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{2}.VAL=1:size(V,1);

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='sx;sy;sz';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
%% Quick viewing of the FEBio input file structure
%The febView function can be used to view the xml structure in a MATLAB figure window.
% 
%disp('Viewing the febio file');
%febView(febio_spec); %Viewing the febio file

%% Exporting the FEBio input file
% Exporting the febio_spec structure to an FEBio input file is done using
% the |febioStruct2xml| function.
febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

% % 
% % 
% % %% Running the FEBio analysis
% % % To run the analysis defined by the created FEBio input file the
% % % |runMonitorFEBio| function is used. The input for this function is a
% % % structure defining job settings e.g. the FEBio input file name. The
% % % optional output runFlag informs the user if the analysis was run
% % % succesfully. 

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; %Display information on the command window
febioAnalysis.disp_log_on=1; %Display convergence information in the command window
febioAnalysis.runMode='internal';%'internal';
febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=3; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!
