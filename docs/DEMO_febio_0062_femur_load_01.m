%% DEMO_febio_0062_femur_load_01
% Below is a demonstration for:
% 
% * Building geometry femur bone
% * Applying a load to the femur head

%% Keywords
%
% * febio_spec version 2.5
% * febio, FEBio
% * beam force loading
% * force control boundary condition
% * tetrahedral elements, tet4
% * femur
% * static, solid
% * displacement logfile
% * stress logfile

%%

clear; close all; clc;

%% Plot settings
fontSize=20;
faceAlpha1=0.8;
markerSize=40;
markerSize2=20;
lineWidth=3;

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');
pathNameSTL=fullfile(defaultFolder,'data','STL'); 
saveName_SED=fullfile(savePath,'SED_no_implant.mat');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stresses
febioLogFileName_strainEnergy=[febioFebFileNamePart,'_energy_out.txt']; %Log file name for exporting strain energy density

%Geometric parameters
distanceCut=250; %Distance from femur to cut bone at
corticalThickness=3; %Thickness used for cortical material definition
volumeFactor=2; %Factor to scale desired volume for interior elements w.r.t. boundary elements

%Define applied force 

forceAbductor=[564.831 -132.696 704.511];
forceVastusLateralis_Walking=[-7.857 -161.505 -811.017];
forceVastusLateralis_StairClimbing=[-19.206 -195.552 -1,179.423];
forceVastusMedialis_StairClimbing=[-76.824 -345.708 -2,331.783]; 
forceVM_inactive=[0 0 0];

n=1;
switch n
    case 1 
        forceVastusLateralis=forceVastusLateralis_Walking; 
        forceVastusMedialis=forceVM_inactive;
    case 2 
        forceVastusLateralis=forceVastusLateralis_StairClimbing;
        forceVastusMedialis=forceVastusMedialis_StairClimbing; 
    otherwise 
        forceVastusLateralis=forceVastusLateralis_StairClimbing;
        forceVastusMedialis=forceVastusMedialis_StairClimbing;
end 
% Distance markers and scaling factor
zLevelWidthMeasure = -75;
zLevelFCML = -395;
scaleFactorSize=1;
distanceMuscleAttachAbductor=15;
distanceMuscleVastusLateralis=10;
distanceMuscleAttachVastusMedialis=10;

%Material parameters (MPa if spatial units are mm)
% Cortical bone
E_youngs1=17000; %Youngs modulus
nu1=0.25; %Poissons ratio

% Cancellous bone
E_youngs2=1500; %Youngs modulus
nu2=0.25; %Poissons ratio

% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
runMode='external'; %'external' or 'internal'

%% Import bone surface model
[stlStruct] = import_STL(fullfile(pathNameSTL,'femur_iso.stl'));
F_bone=stlStruct.solidFaces{1}; %Faces
V_bone=stlStruct.solidVertices{1}; %Vertices

%% Scale and reorient

V_bone=V_bone.*1000; %Scale to mm
V_bone=V_bone.*scaleFactorSize; %Scale size further

[F_bone,V_bone]=mergeVertices(F_bone,V_bone); % Merging nodes
Q1=euler2DCM([0 0 0.065*pi]);
V_bone=V_bone*Q1;
Q2=euler2DCM([-0.5*pi 0 0]);
V_bone=V_bone*Q2;
Q3=euler2DCM([0 0 0.36*pi]);
V_bone=V_bone*Q3;

%% Visualize bone surface

cFigure; hold on;
gpatch(F_bone,V_bone,'w','k',1);
% patchNormPlot(F_bone,V_bone)
axisGeom; camlight headlight;
drawnow;

%% Determining femur length and person's height 

% The Estimation of Stature on the Basis of Measurements of the Femur - height 
% Estimating body mass and composition from proximal femur dimensions using dual energy x-ray absorptiometry - weight 

% P3 = [-15 20 -75];
% [~,indNode3]=minDist(P3,V_bone);

femurLength = max(V_bone(:,3)) - min(V_bone(:,3)); %femur length 

[~,V_bone_slice,~,~,Eb]=triSurfSlice(F_bone,V_bone,[],[0 0 scaleFactorSize.*zLevelWidthMeasure],[0 0 1]);
indSliceCurve=edgeListToCurve(Eb);
V_slice_curve=V_bone_slice(indSliceCurve,:);
sliceArea = polyarea(V_slice_curve(:,1),V_slice_curve(:,2));
subtrochanterMedLatDia = sqrt(sliceArea./(0.25*pi));

[~,V_bone_slice_bottom,~,~,Eb2]=triSurfSlice(F_bone,V_bone,[],[0 0 scaleFactorSize.*zLevelFCML],[0 0 1]);
indSliceCurveBottom=edgeListToCurve(Eb2);
V_slice_curve_bottom=V_bone_slice_bottom(indSliceCurveBottom,:);
sliceAreaBottom = polyarea(V_slice_curve_bottom(:,1),V_slice_curve_bottom(:,1));

FCML = max(V_bone_slice_bottom(:,1)) - min(V_bone_slice_bottom(:,1)); %mediolateral breadth of the articular surface of the femoral condyles

% subtrochanterMedLatDia = distND(V_bone(indNode3,:),V_bone(indNode4,:)); %subtrochanter medio-latreal diameter

bodyHeight = 4*femurLength;
bodyMass = 1.37 * (FCML-42.8);

forceTotal=[-0.54 -0.32  -2.292].*bodyMass;

%%

cFigure; hold on;
gpatch(F_bone,V_bone,'w','none',0.5);
plotV(V_bone_slice(indSliceCurve,:),'r-','LineWidth',5)
plotV(V_bone_slice_bottom(indSliceCurveBottom,:),'r-','LineWidth',5)
axisGeom; 
camlight headlight;
drawnow;

%% Cut bone surface

%Slicing surface

[F_bone,V_bone,~,logicSide,~]=triSurfSlice(F_bone,V_bone,[],[0 0 -distanceCut],[0 0 1]);

F_bone=F_bone(logicSide==0,:);
[F_bone,V_bone]=patchCleanUnused(F_bone,V_bone);

Eb=patchBoundary(F_bone,V_bone);
indCurve=edgeListToCurve(Eb);
indCurve=indCurve(1:end-1);

cparSmooth.n=5;
cparSmooth.Method='HC';
[V_Eb_smooth]=patchSmooth(Eb,V_bone(:,[1 2]),[],cparSmooth);
V_bone(indCurve,[1 2])=V_Eb_smooth(indCurve,:);

cparSmooth.n=5;
cparSmooth.Method='HC';
cparSmooth.RigidConstraints=indCurve;
[V_bone]=patchSmooth(F_bone,V_bone,[],cparSmooth);

pointSpacing=mean(patchEdgeLengths(F_bone,V_bone));

[F_bone2,V_bone2]=regionTriMesh3D({V_bone(indCurve,:)},pointSpacing,0,'linear');
if dot(mean(patchNormal(F_bone2,V_bone2)),[0 0 1])>0
    F_bone2=fliplr(F_bone2);
end

[F_bone,V_bone,C_bone]=joinElementSets({F_bone,F_bone2},{V_bone,V_bone2});
[F_bone,V_bone]=mergeVertices(F_bone,V_bone);

%% Visualize bone surface

cFigure; hold on;
gpatch(F_bone,V_bone,C_bone,'k',1);
patchNormPlot(F_bone,V_bone);
axisGeom; camlight headlight;
drawnow;

%% Mesh using tetgen

%Find interior point
V_inner_bone=getInnerPoint(F_bone,V_bone);

%% Visualize interior point
cFigure; hold on;
gpatch(F_bone,V_bone,'w','none',0.5);
plotV(V_inner_bone,'r.','MarkerSize',25)
axisGeom; camlight headlight;
drawnow;

%% 
% Regional mesh volume parameter
tetVolume=tetVolMeanEst(F_bone,V_bone); %Volume for regular tets

tetGenStruct.stringOpt='-pq1.2AaY';
tetGenStruct.Faces=F_bone;
tetGenStruct.Nodes=V_bone;
tetGenStruct.holePoints=[];
tetGenStruct.faceBoundaryMarker=C_bone; %Face boundary markers
tetGenStruct.regionPoints=V_inner_bone; %region points
tetGenStruct.regionA=tetVolume*volumeFactor;

[meshOutput]=runTetGen(tetGenStruct); %Run tetGen 

% Access elements, nodes, and boundary faces
E=meshOutput.elements;
V=meshOutput.nodes;
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
CE=meshOutput.elementMaterialID;

%% Define material regions in bone

indBoundary=unique(Fb(Cb==1,:));
DE=minDist(V,V(indBoundary,:));
logicCorticalNodes=DE<=corticalThickness; 
logicCorticalElements=any(logicCorticalNodes(E),2);
logicCancellousElements=~logicCorticalElements;

E1=E(logicCorticalElements,:);
E2=E(logicCancellousElements,:);
E=[E1;E2];
elementMaterialID=[ones(size(E1,1),1);2*ones(size(E2,1),1);];
meshOutput.elements=E;
meshOutput.elementMaterialID=elementMaterialID;

%% Visualizing solid mesh 

hFig=cFigure; hold on;
optionStruct.hFig=hFig;
meshView(meshOutput,optionStruct);
axisGeom; 
drawnow;

%% Find femoral head

w=100;
f=[1 2 3 4];
v=w*[-1 -1 0; -1 1 0; 1 1 0; 1 -1 0];

p=[0 0 0];
Q=euler2DCM([0 (150/180)*pi 0]);
v=v*Q;
v=v+p;

Vr=V*Q';
Vr=Vr+p;
logicHeadNodes=Vr(:,3)<0;
logicHeadFaces=all(logicHeadNodes(Fb),2);
bcPrescribeList=unique(Fb(logicHeadFaces,:));

%%
% Visualize femoral head nodes for prescribed force boundary conditions
cFigure; 
hold on;
gpatch(Fb,V,'w','k',1);
gpatch(f,v,'r','k',0.5);
plotV(V(bcPrescribeList,:),'r.','markerSize',15)
axisGeom; camlight headlight;
drawnow;

%% Work out force distribution on femoral head surface nodes
% This is based on surface normal directions. Forces are assumed to only be
% able to act in a compressive sense on the bone. 

[~,~,N]=patchNormal(fliplr(Fb),V); %Nodal normal directions

FX=[forceTotal(1) 0 0]; %X force vector
FY=[0 forceTotal(2) 0]; %Y force vector
FZ=[0 0 forceTotal(3)]; %Z force vector

wx=dot(N(bcPrescribeList,:),FX(ones(numel(bcPrescribeList),1),:),2);
wy=dot(N(bcPrescribeList,:),FY(ones(numel(bcPrescribeList),1),:),2);
wz=dot(N(bcPrescribeList,:),FZ(ones(numel(bcPrescribeList),1),:),2);

%Force zero
wx(wx>0)=0; wy(wy>0)=0; wz(wz>0)=0;

force_X=forceTotal(1).*ones(numel(bcPrescribeList),1).*wx;
force_Y=forceTotal(2).*ones(numel(bcPrescribeList),1).*wy;
force_Z=forceTotal(3).*ones(numel(bcPrescribeList),1).*wz;

force_X=force_X./sum(force_X(:)); %sum now equal to 1
force_X=force_X.*forceTotal(1); %sum now equal to desired

force_Y=force_Y./sum(force_Y(:)); %sum now equal to 1
force_Y=force_Y.*forceTotal(2); %sum now equal to desired

force_Z=force_Z./sum(force_Z(:)); %sum now equal to 1
force_Z=force_Z.*forceTotal(3); %sum now equal to desired

%%
cFigure; 
subplot(1,3,1);hold on;
title('F_x');
gpatch(Fb,V,'w','none',0.5);
quiverVec([0 0 0],FX,100,'k');
% scatterV(V(indicesHeadNodes,:),15)
quiverVec(V(bcPrescribeList,:),N(bcPrescribeList,:),10,force_X);
axisGeom; camlight headlight;
colormap(gca,gjet(250)); colorbar; 

subplot(1,3,2);hold on;
title('F_y');
gpatch(Fb,V,'w','none',0.5);
quiverVec([0 0 0],FY,100,'k');
% scatterV(V(indicesHeadNodes,:),15)
quiverVec(V(bcPrescribeList,:),N(bcPrescribeList,:),10,force_Y);
axisGeom; camlight headlight;
colormap(gca,gjet(250)); colorbar; 

subplot(1,3,3);hold on;
title('F_z');
gpatch(Fb,V,'w','none',0.5);
quiverVec([0 0 0],FZ,100,'k');
% scatterV(V(indicesHeadNodes,:),15)
quiverVec(V(bcPrescribeList,:),N(bcPrescribeList,:),10,force_Z);
axisGeom; camlight headlight;
colormap(gca,gjet(250)); colorbar; 

drawnow;

%% Marking muscle locations

P_abductor_find = [-69.771045288206111 8.185179717034659 -5.575329878303917]; %Coordinate at centre of muscle attachment
[~,indAbductor]=minDist(P_abductor_find,V); %Node number of point at centre of attachment
dAbductor=meshDistMarch(Fb,V,indAbductor); %Distance (on mesh) from attachement centre
bcPrescibeList_abductor=find(dAbductor<=distanceMuscleAttachAbductor); %Node numbers for attachment site

P_VastusLateralis_find = [-58.763839901506827 19.145444610053566 -51.005278396808819]; %Coordinate at centre of muscle attachment
[~,indVastusLateralis]=minDist(P_VastusLateralis_find,V); %Node number of point at centre of attachment
dVastusLateralis=meshDistMarch(Fb,V,indVastusLateralis); %Distance (on mesh) from attachement centre
bcPrescibeList_VastusLateralis=find(dVastusLateralis<=distanceMuscleVastusLateralis); %Node numbers for attachment site


P_VastusMedialis_find = [-18.533631492778085 9.501312355952791 -85.666499329588035]; %Coordinate at centre of muscle attachment
[~,indVastusMedialis]=minDist(P_VastusMedialis_find,V); %Node number of point at centre of attachment
dVastusMedialis=meshDistMarch(Fb,V,indVastusMedialis); %Distance (on mesh) from attachement centre
bcPrescibeList_VastusMedialis=find(dVastusMedialis<=distanceMuscleAttachVastusMedialis); %Node numbers for attachment site

%%  Muscle force definition 
forceAbductor_distributed=forceAbductor.*ones(numel(bcPrescibeList_abductor),1)./numel(bcPrescibeList_abductor);

forceVastusLateralis_distributed=forceVastusLateralis.*ones(numel(bcPrescibeList_VastusLateralis),1)./numel(bcPrescibeList_VastusLateralis);

forceVastusMedialis_distributed=forceVastusMedialis.*ones(numel(bcPrescibeList_VastusMedialis),1)./numel(bcPrescibeList_VastusMedialis);

%% Visualizing boundary conditions

F_bottomSupport=Fb(Cb==2,:);
bcSupportList=unique(F_bottomSupport(:));

cFigure; hold on;
gpatch(Fb,V,'kw','none',0.7);
hl(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',25);
hl(2)=plotV(V(bcPrescribeList,:),'r.','MarkerSize',25);
hl(3)=plotV(V(bcPrescibeList_abductor,:),'g.','MarkerSize',25);
hl(4)=plotV(V(bcPrescibeList_VastusLateralis,:),'b.','MarkerSize',25);
hl(5)=plotV(V(bcPrescibeList_VastusMedialis,:),'g.','MarkerSize',25);
legend(hl,{'BC support','BC force prescribe','MAP abductor','MAP Vastus Lateralis','MAP Vastus Medialis'});
axisGeom; 
camlight headlight;
drawnow;

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='2.5'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

%Control section
febio_spec.Control.analysis.ATTR.type='static';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;
febio_spec.Control.max_refs=max_refs;
febio_spec.Control.max_ups=max_ups;

%Material section
febio_spec.Material.material{1}.ATTR.type='neo-Hookean';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.E=E_youngs1;
febio_spec.Material.material{1}.v=nu1;

febio_spec.Material.material{2}.ATTR.type='neo-Hookean';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.E=E_youngs2;
febio_spec.Material.material{2}.v=nu2;

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='tet4'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name='CorticalBone'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(E1,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=E1;

febio_spec.Geometry.Elements{2}.ATTR.type='tet4'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat=2; %material index for this set 
febio_spec.Geometry.Elements{2}.ATTR.name='CancellousBone'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL=E2;

% -> NodeSets
febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcSupportList(:);

febio_spec.Geometry.NodeSet{2}.ATTR.name='indicesHeadSurfaceNodes';
febio_spec.Geometry.NodeSet{2}.node.ATTR.id=bcPrescribeList(:);

febio_spec.Geometry.NodeSet{3}.ATTR.name='indicesAbductor';
febio_spec.Geometry.NodeSet{3}.node.ATTR.id=bcPrescibeList_abductor(:);

febio_spec.Geometry.NodeSet{4}.ATTR.name='indicesVastusLateralis';
febio_spec.Geometry.NodeSet{4}.node.ATTR.id=bcPrescibeList_VastusLateralis(:);

febio_spec.Geometry.NodeSet{5}.ATTR.name='indicesVastusMedialis';
febio_spec.Geometry.NodeSet{5}.node.ATTR.id=bcPrescibeList_VastusMedialis(:);

%Boundary condition section
% -> Fix boundary conditions
febio_spec.Boundary.fix{1}.ATTR.bc='x';
febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{2}.ATTR.bc='y';
febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{3}.ATTR.bc='z';
febio_spec.Boundary.fix{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;

febio_spec.MeshData.NodeData{1}.ATTR.name='force_X';
febio_spec.MeshData.NodeData{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.MeshData.NodeData{1}.node.VAL=force_X;
febio_spec.MeshData.NodeData{1}.node.ATTR.lid=(1:1:numel(bcPrescribeList))';

febio_spec.MeshData.NodeData{2}.ATTR.name='force_Y';
febio_spec.MeshData.NodeData{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.MeshData.NodeData{2}.node.VAL=force_Y;
febio_spec.MeshData.NodeData{2}.node.ATTR.lid=(1:1:numel(bcPrescribeList))';

febio_spec.MeshData.NodeData{3}.ATTR.name='force_Z';
febio_spec.MeshData.NodeData{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.MeshData.NodeData{3}.node.VAL=force_Z;
febio_spec.MeshData.NodeData{3}.node.ATTR.lid=(1:1:numel(bcPrescribeList))';

febio_spec.MeshData.NodeData{4}.ATTR.name='forceAbductor_X';
febio_spec.MeshData.NodeData{4}.ATTR.node_set=febio_spec.Geometry.NodeSet{3}.ATTR.name;
febio_spec.MeshData.NodeData{4}.node.VAL=forceAbductor_distributed(:,1);
febio_spec.MeshData.NodeData{4}.node.ATTR.lid=(1:1:numel(bcPrescibeList_abductor))';

febio_spec.MeshData.NodeData{5}.ATTR.name='forceAbductor_Y';
febio_spec.MeshData.NodeData{5}.ATTR.node_set=febio_spec.Geometry.NodeSet{3}.ATTR.name;
febio_spec.MeshData.NodeData{5}.node.VAL=forceAbductor_distributed(:,2);
febio_spec.MeshData.NodeData{5}.node.ATTR.lid=(1:1:numel(bcPrescibeList_abductor))';

febio_spec.MeshData.NodeData{6}.ATTR.name='forceAbductor_Z';
febio_spec.MeshData.NodeData{6}.ATTR.node_set=febio_spec.Geometry.NodeSet{3}.ATTR.name;
febio_spec.MeshData.NodeData{6}.node.VAL=forceAbductor_distributed(:,3);
febio_spec.MeshData.NodeData{6}.node.ATTR.lid=(1:1:numel(bcPrescibeList_abductor))';

febio_spec.MeshData.NodeData{7}.ATTR.name='forceVL_X';
febio_spec.MeshData.NodeData{7}.ATTR.node_set=febio_spec.Geometry.NodeSet{4}.ATTR.name;
febio_spec.MeshData.NodeData{7}.node.VAL=forceVastusLateralis_distributed(:,1);
febio_spec.MeshData.NodeData{7}.node.ATTR.lid=(1:1:numel(bcPrescibeList_VastusLateralis))';

febio_spec.MeshData.NodeData{8}.ATTR.name='forceVL_Y';
febio_spec.MeshData.NodeData{8}.ATTR.node_set=febio_spec.Geometry.NodeSet{4}.ATTR.name;
febio_spec.MeshData.NodeData{8}.node.VAL=forceVastusLateralis_distributed(:,2);
febio_spec.MeshData.NodeData{8}.node.ATTR.lid=(1:1:numel(bcPrescibeList_VastusLateralis))';

febio_spec.MeshData.NodeData{9}.ATTR.name='forceVL_Z';
febio_spec.MeshData.NodeData{9}.ATTR.node_set=febio_spec.Geometry.NodeSet{4}.ATTR.name;
febio_spec.MeshData.NodeData{9}.node.VAL=forceVastusLateralis_distributed(:,3);
febio_spec.MeshData.NodeData{9}.node.ATTR.lid=(1:1:numel(bcPrescibeList_VastusLateralis))';

febio_spec.MeshData.NodeData{10}.ATTR.name='forceVM_X';
febio_spec.MeshData.NodeData{10}.ATTR.node_set=febio_spec.Geometry.NodeSet{5}.ATTR.name;
febio_spec.MeshData.NodeData{10}.node.VAL=forceVastusMedialis_distributed(:,1);
febio_spec.MeshData.NodeData{10}.node.ATTR.lid=(1:1:numel(bcPrescibeList_VastusMedialis))';

febio_spec.MeshData.NodeData{11}.ATTR.name='forceVM_Y';
febio_spec.MeshData.NodeData{11}.ATTR.node_set=febio_spec.Geometry.NodeSet{5}.ATTR.name;
febio_spec.MeshData.NodeData{11}.node.VAL=forceVastusMedialis_distributed(:,2);
febio_spec.MeshData.NodeData{11}.node.ATTR.lid=(1:1:numel(bcPrescibeList_VastusMedialis))';

febio_spec.MeshData.NodeData{12}.ATTR.name='forceVM_Z';
febio_spec.MeshData.NodeData{12}.ATTR.node_set=febio_spec.Geometry.NodeSet{5}.ATTR.name;
febio_spec.MeshData.NodeData{12}.node.VAL=forceVastusMedialis_distributed(:,3);
febio_spec.MeshData.NodeData{12}.node.ATTR.lid=(1:1:numel(bcPrescibeList_VastusMedialis))';

%Loads section
% -> Prescribed nodal forces
febio_spec.Loads.nodal_load{1}.ATTR.bc='x';
febio_spec.Loads.nodal_load{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Loads.nodal_load{1}.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load{1}.scale.VAL=1;
febio_spec.Loads.nodal_load{1}.value.ATTR.node_data=febio_spec.MeshData.NodeData{1}.ATTR.name;

febio_spec.Loads.nodal_load{2}.ATTR.bc='y';
febio_spec.Loads.nodal_load{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Loads.nodal_load{2}.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load{2}.scale.VAL=1;
febio_spec.Loads.nodal_load{2}.value.ATTR.node_data=febio_spec.MeshData.NodeData{2}.ATTR.name;

febio_spec.Loads.nodal_load{3}.ATTR.bc='z';
febio_spec.Loads.nodal_load{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Loads.nodal_load{3}.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load{3}.scale.VAL=1;
febio_spec.Loads.nodal_load{3}.value.ATTR.node_data=febio_spec.MeshData.NodeData{3}.ATTR.name;

febio_spec.Loads.nodal_load{4}.ATTR.bc='x';
febio_spec.Loads.nodal_load{4}.ATTR.node_set=febio_spec.Geometry.NodeSet{3}.ATTR.name;
febio_spec.Loads.nodal_load{4}.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load{4}.scale.VAL=1;
febio_spec.Loads.nodal_load{4}.value.ATTR.node_data=febio_spec.MeshData.NodeData{4}.ATTR.name;

febio_spec.Loads.nodal_load{5}.ATTR.bc='y';
febio_spec.Loads.nodal_load{5}.ATTR.node_set=febio_spec.Geometry.NodeSet{3}.ATTR.name;
febio_spec.Loads.nodal_load{5}.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load{5}.scale.VAL=1;
febio_spec.Loads.nodal_load{5}.value.ATTR.node_data=febio_spec.MeshData.NodeData{5}.ATTR.name;

febio_spec.Loads.nodal_load{6}.ATTR.bc='z';
febio_spec.Loads.nodal_load{6}.ATTR.node_set=febio_spec.Geometry.NodeSet{3}.ATTR.name;
febio_spec.Loads.nodal_load{6}.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load{6}.scale.VAL=1;
febio_spec.Loads.nodal_load{6}.value.ATTR.node_data=febio_spec.MeshData.NodeData{6}.ATTR.name;

febio_spec.Loads.nodal_load{7}.ATTR.bc='x';
febio_spec.Loads.nodal_load{7}.ATTR.node_set=febio_spec.Geometry.NodeSet{4}.ATTR.name;
febio_spec.Loads.nodal_load{7}.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load{7}.scale.VAL=1;
febio_spec.Loads.nodal_load{7}.value.ATTR.node_data=febio_spec.MeshData.NodeData{7}.ATTR.name;

febio_spec.Loads.nodal_load{8}.ATTR.bc='y';
febio_spec.Loads.nodal_load{8}.ATTR.node_set=febio_spec.Geometry.NodeSet{4}.ATTR.name;
febio_spec.Loads.nodal_load{8}.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load{8}.scale.VAL=1;
febio_spec.Loads.nodal_load{8}.value.ATTR.node_data=febio_spec.MeshData.NodeData{8}.ATTR.name;

febio_spec.Loads.nodal_load{9}.ATTR.bc='z';
febio_spec.Loads.nodal_load{9}.ATTR.node_set=febio_spec.Geometry.NodeSet{4}.ATTR.name;
febio_spec.Loads.nodal_load{9}.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load{9}.scale.VAL=1;
febio_spec.Loads.nodal_load{9}.value.ATTR.node_data=febio_spec.MeshData.NodeData{9}.ATTR.name;

febio_spec.Loads.nodal_load{10}.ATTR.bc='x';
febio_spec.Loads.nodal_load{10}.ATTR.node_set=febio_spec.Geometry.NodeSet{5}.ATTR.name;
febio_spec.Loads.nodal_load{10}.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load{10}.scale.VAL=1;
febio_spec.Loads.nodal_load{10}.value.ATTR.node_data=febio_spec.MeshData.NodeData{10}.ATTR.name;

febio_spec.Loads.nodal_load{11}.ATTR.bc='y';
febio_spec.Loads.nodal_load{11}.ATTR.node_set=febio_spec.Geometry.NodeSet{5}.ATTR.name;
febio_spec.Loads.nodal_load{11}.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load{11}.scale.VAL=1;
febio_spec.Loads.nodal_load{11}.value.ATTR.node_data=febio_spec.MeshData.NodeData{11}.ATTR.name;

febio_spec.Loads.nodal_load{12}.ATTR.bc='z';
febio_spec.Loads.nodal_load{12}.ATTR.node_set=febio_spec.Geometry.NodeSet{5}.ATTR.name;
febio_spec.Loads.nodal_load{12}.scale.ATTR.lc=1;
febio_spec.Loads.nodal_load{12}.scale.VAL=1;
febio_spec.Loads.nodal_load{12}.value.ATTR.node_data=febio_spec.MeshData.NodeData{12}.ATTR.name;

%Output section 
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='s1;s2;s3';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.element_data{1}.VAL=1:1:size(E,1); %Rigid body material id

febio_spec.Output.logfile.element_data{2}.ATTR.file=febioLogFileName_strainEnergy;
febio_spec.Output.logfile.element_data{2}.ATTR.data='sed';
febio_spec.Output.logfile.element_data{2}.ATTR.delim=',';
febio_spec.Output.logfile.element_data{2}.VAL=1:1:size(E,1);

%% Quick viewing of the FEBio input file structure
% The |febView| function can be used to view the xml structure in a MATLAB
% figure window. 

%%
% |febView(febio_spec); %Viewing the febio file|

%% Exporting the FEBio input file
% Exporting the febio_spec structure to an FEBio input file is done using
% the |febioStruct2xml| function. 

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

%% Running the FEBio analysis
% To run the analysis defined by the created FEBio input file the
% |runMonitorFEBio| function is used. The input for this function is a
% structure defining job settings e.g. the FEBio input file name. The
% optional output runFlag informs the user if the analysis was run
% succesfully. 

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; %Display information on the command window
febioAnalysis.disp_log_on=1; %Display convergence information in the command window
febioAnalysis.runMode=runMode;%'internal';
febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results 

if runFlag==1 %i.e. a succesful run
    
    % Importing nodal displacement from a log file
    [time_mat, N_disp_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp)); %Nodal displacement    
    time_mat=[0; time_mat(:)]; %Time

    N_disp_mat=N_disp_mat(:,2:end,:);
    sizImport=size(N_disp_mat);
    sizImport(3)=sizImport(3)+1;
    N_disp_mat_n=zeros(sizImport);
    N_disp_mat_n(:,:,2:end)=N_disp_mat;
    N_disp_mat=N_disp_mat_n;
    DN=N_disp_mat(:,:,end);
    DN_magnitude=sqrt(sum(DN(:,3).^2,2));
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);    

    %%
    % Importing element strain energies from a log file
    [~,E_energy,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_strainEnergy)); %Element strain energy
    
    %Remove nodal index column
    E_energy=E_energy(:,2:end,:);
    
    %Add initial state i.e. zero energy
    sizImport=size(E_energy);
    sizImport(3)=sizImport(3)+1;
    E_energy_mat_n=zeros(sizImport);
    E_energy_mat_n(:,:,2:end)=E_energy;
    E_energy=E_energy_mat_n;
   
    %%
    [FE_face,C_energy_face]=element2patch(E,E_energy(:,:,end),'tet4');
    [CV]=faceToVertexMeasure(FE_face,V,C_energy_face);
    [indBoundary]=tesBoundary(FE_face,V);
    Fb=FE_face(indBoundary,:);

    %% Saving strain energy data for comparison in implant demo
    
    %Element centre coordinates
    VE=patchCentre(E,V);
    
    %Construct interpolation function
    interpFuncEnergy=scatteredInterpolant(VE,E_energy(:,:,end),'natural','linear');
    
    %Store components in structure
    outputStruct.SED=E_energy;
    outputStruct.elements=E;
    outputStruct.nodes=V;
    outputStruct.boundaryFaces=Fb;
    outputStruct.elementCenters=VE;
    outputStruct.interpFuncEnergy=interpFuncEnergy;
    
    %Save structure
    save(saveName_SED,'-struct','outputStruct');
            
    %%
    
    C_data=E_energy(:,:,end);
    
    n=[0 1 0]; %Normal direction to plane
    P=mean(V,1); %Point on plane
    [logicAt]=meshCleave(E,V,P,n);
    
    % Get faces and matching color data for visualization
    [F_cleave,CF_cleave]=element2patch(E(logicAt,:),C_data(logicAt));
    
    hf=cFigure; hold on;
    title('Strain energy density');
    gpatch(Fb,V,'w','none',0.1); %Add graphics object to animate
    
    hp1=gpatch(F_cleave,V,CF_cleave,'k',1);
    axisGeom(gca,fontSize);
    colormap(gjet(250));
    caxis([0 max(C_data)]/2);
    colorbar; axis manual;
    camlight headligth;
    gdrawnow;
    
    nSteps=25; %Number of animation steps
    
    %Create the time vector
    animStruct.Time=linspace(0,1,nSteps);
    
    %The vector lengths
    y=linspace(min(V(:,2)),max(V(:,2)),nSteps);
    for q=1:1:nSteps
        %Get logic for slice of elements
        logicAt=meshCleave(E,V,[0 y(q) 0],n,[1 0]);
        
        % Get faces and matching color data for cleaves elements
        [F_cleave,CF_cleave]=element2patch(E(logicAt,:),C_data(logicAt));
        
        %Set entries in animation structure
        animStruct.Handles{q}=[hp1 hp1]; %Handles of objects to animate
        animStruct.Props{q}={'Faces','CData'}; %Properties of objects to animate
        animStruct.Set{q}={F_cleave,CF_cleave}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct);
    
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    axLim=[min(min(V_DEF,[],3),[],1); max(max(V_DEF,[],3),[],1)];
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure  
    title('Strain energy density')
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1=gpatch(Fb,V_DEF(:,:,end),CV,'k',1); %Add graphics object to animate
    hp1.FaceColor='Interp';
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([0 max(E_energy(:))/25]);    
    axis(axLim(:)'); %Set axis limits statically    
    camlight headlight;        
    
    % Set up animation features
    animStruct.Time=time_mat; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        DN=N_disp_mat(:,:,qt); %Current displacancellousBone

        [FE_face,C_energy_face]=element2patch(E,E_energy(:,:,qt),'tet4');
        [CV]=faceToVertexMeasure(FE_face,V,C_energy_face);
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1 hp1]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;

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
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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