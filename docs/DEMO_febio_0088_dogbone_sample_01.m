%% DEMO_febio_0001_cube_uniaxial
% Below is a demonstration for:
% 
% * Building geometry for a cube with hexahedral elements
% * Defining the boundary conditions 
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement and stress results

%% Keywords
%
% * febio_spec version 4.0
% * febio, FEBio
% * uniaxial loading
% * compression, tension, compressive, tensile
% * displacement control, displacement boundary condition
% * hexahedral elements, hex8
% * cube, box, rectangular
% * static, solid
% * neo-hookean, uncoupled
% * displacement logfile
% * stress logfile

%%

clear; close all; clc;

%% Plot settings
fontSize=20;
faceAlpha1=0.8;
markerSize=40;
markerSize2=35;
lineWidth=3;
cMap=viridis(20); %colormap 

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath'))); 
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_strain=[febioFebFileNamePart,'_energy_out.txt']; %Log file name for exporting strain energy density

%Specifying dimensions and number of elements
centreWidth = 6;
centreLength = 33;
fullLength = 115;
tabAdd = 9.5;
alphaTab = pi/4; 
r_tab = 9;
sampleThickness = 2.25; 
tabWidth = centreWidth+2*tabAdd;
pointSpacing = 1.5; 
pointSpacingCentre = pointSpacing/2; 
pointSpacingThickness = pointSpacing/2; 
sampleGripGripHeight = 68; 

%Loading parameters
appliedLinearStrain=0.3;
clampCompressiveLinearStrain=0.2;

%Initial material parameter set
c1=1e-3;
m1=2;
k_factor=100;
k=c1*k_factor;

% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

runMode='internal';% 'internal' or 'external'

%%

clampCompressiveDisplacement=(sampleThickness.*clampCompressiveLinearStrain)/2;
clampTensionDisplacement=(sampleGripGripHeight.*appliedLinearStrain);

%% Create quarter segment boundary curve

Vc1 = [ 0 sampleGripGripHeight/2 0; ...
        0 fullLength/2 0; ...
       -tabWidth/2 fullLength/2 0; % top left
       -tabWidth/2 sampleGripGripHeight/2 0;           
      ];
Vc1=flipud(Vc1);

optionStruct.w = sampleGripGripHeight/2-centreLength/2;
optionStruct.h = tabAdd;
optionStruct.p = 1; 
optionStruct.n = 100; 
[Vc2]=roundedSigmoid(optionStruct);
Vc2=Vc2(:,[2,1]);
Vc2(:,2)=-Vc2(:,2);
Vc2(:,1)=Vc2(:,1)+optionStruct.h/2-tabWidth/2;
Vc2(:,2)=Vc2(:,2)-optionStruct.w/2+sampleGripGripHeight/2;
Vc2(:,3)=0;

Vc3 =[Vc2(end,:);...
    -centreWidth/2,0.0,0.0;...
      0.0,0.0,0.0;...
      0.0 Vc2(end,2) 0];

Vc1 = evenlySpaceCurve(Vc1,pointSpacing,'linear',0,1:size(Vc1,1));
Vc2 = evenlySpaceCurve(Vc2,pointSpacing,'linear',0);
Vc3 = evenlySpaceCurve(Vc3,pointSpacingCentre,'linear',0,1:size(Vc3,1));
Vc4 = evenlySpaceCurve([0.0 sampleGripGripHeight/2 0; -tabWidth/2 sampleGripGripHeight/2 0],pointSpacing,'linear',0);
Vc5 = evenlySpaceCurve([0.0 Vc2(end,2) 0.0; 0.0 sampleGripGripHeight/2 0],pointSpacing,'linear',0);
Vc6 = evenlySpaceCurve([Vc2(end,:); 0.0 Vc2(end,2) 0.0],pointSpacingCentre,'linear',0);
[Fs11,Vs11] = regionTriMesh2D({[Vc4(2:end-1,[1,2]); Vc1(:,[1,2])]},pointSpacing,0,0);
[Fs12,Vs12] = regionTriMesh2D({[Vc2(:,[1,2]); Vc6(2:end,[1,2]);Vc5(2:end,[1,2]);Vc4(2:end-1,[1,2])]},pointSpacing,0,0);
[Fs13,Vs13] = regionTriMesh2D({[Vc3(:,[1,2]); flipud(Vc6(2:end-1,[1,2]));]},pointSpacingCentre,0,0);
[Fs,Vs,Cs] = joinElementSets({Fs11,Fs12,Fs13},{Vs11,Vs12,Vs13});
[Fs,Vs] = mergeVertices(Fs,Vs);
Vs(:,3)=0.0;

cFigure; hold on; 
gpatch(Fs,Vs,Cs,'k');
% gpatch(Fs11,Vs11,'rw','k');
% gpatch(Fs12,Vs12,'gw','k');
% gpatch(Fs13,Vs13,'bw','k');

plotV(Vc1,'b.-','MarkerSize',15,'LineWidth',2)
plotV(Vc2,'r.-','MarkerSize',15,'LineWidth',2)
plotV(Vc3,'g.-','MarkerSize',15,'LineWidth',2)
plotV(Vc4,'y.-','MarkerSize',15,'LineWidth',2)
plotV(Vc5,'c.-','MarkerSize',15,'LineWidth',2)
plotV(Vc6,'k.-','MarkerSize',15,'LineWidth',2)
plotV(Vc1(1,:),'b.-','MarkerSize',50,'LineWidth',2)
plotV(Vc2(1,:),'r.-','MarkerSize',50,'LineWidth',2)

axisGeom(gca,fontSize); camlight headlight; 
gdrawnow;

%%

Fs2 = fliplr(Fs);
Vs2 = Vs;
Vs2(:,1) = -Vs2(:,1);

Fs3 = fliplr(Fs);
Vs3 = Vs;
Vs3(:,2) = -Vs3(:,2);

Fs4 = Fs;
Vs4 = Vs;
Vs4(:,1) = -Vs4(:,1);
Vs4(:,2) = -Vs4(:,2);

[Fs,Vs,Cs] = joinElementSets({Fs,Fs2,Fs3,Fs4},{Vs,Vs2,Vs3,Vs4},{Cs,Cs,Cs,Cs});
[Fs,Vs] = mergeVertices(Fs,Vs);
Vs(:,3) = -sampleThickness/2;

nThick = ceil(sampleThickness/pointSpacingThickness);
[E,V,Fp1,Fp2] =patchThick(Fs,Vs,[0 0 1],sampleThickness,nThick);
Q = euler2DCM([0.5*pi,0,0]); 
V = V*Q; % Rotate upward
CE = repmat(Cs,nThick,1);
[F,CF,CF_type] = element2patch(E,CE,'penta6');
indB = tesBoundary(F);
Fb = {F{1}(indB{1},:),F{2}(indB{2},:)};
Cb = {CF{1}(indB{1},:),CF{2}(indB{2},:)};

logicMiddleElements = CE ==3;

%% 
% Plotting model boundary surfaces and a cut view

hFig=cFigure; hold on; 
title('Thickened solid mesh','FontSize',fontSize);

gpatch(Fb,V,Cb,'k',1); 

axisGeom(gca,fontSize); camlight headlight; 
colormap(gjet(3)); icolorbar; 
gdrawnow;

%% Define clamping surfaces

ny =[0,1,0];
Ft = Fb{1}; % Boundary triangles (top and bottom faces)
Ct = Cb{1};
Nt = patchNormal(Ft,V);
d = dot(Nt,ny(ones(size(Nt,1),1),:),2);
ZF = patchCentre(Ft,V(:,3));

Fc1 = Ft(d<0 & ZF>0 & Ct==1,:);
Fc2 = Ft(d>0 & ZF>0 & Ct==1,:);
Fc3 = Ft(d>0 & ZF<0 & Ct==1,:);
Fc4 = Ft(d<0 & ZF<0 & Ct==1,:);

%% 
% Visualize clamping surfaces 

cFigure; hold on;
title('Clamping surfaces','FontSize',fontSize);
gpatch(Fb,V,'w','none',0.25);
% patchNormPlot(Ft,V)
hp(1)=gpatch(Fc1,V,'rw','k',1);
hp(2)=gpatch(Fc2,V,'gw','k',1);
hp(3)=gpatch(Fc3,V,'bw','k',1);
hp(4)=gpatch(Fc4,V,'yw','k',1);
legend(hp,{'Surf. 1','Surf. 2','Surf. 3','Surf. 4'})
axisGeom(gca,fontSize);
camlight headlight;
drawnow;  

%% Define BC's

bcPrescribeList1=unique(Fc1(:)); % Nodes of surface 1
bcPrescribeList2=unique(Fc2(:)); % Nodes of surface 2
bcPrescribeList3=unique(Fc3(:)); % Nodes of surface 3
bcPrescribeList4=unique(Fc4(:)); % Nodes of surface 4

%%
% Visualize boundary conditions

cFigure; hold on;
title('Complete model','FontSize',fontSize);

gpatch(Fb,V,'w','none',0.25);

hp(1)=plotV(V(bcPrescribeList1,:),'r.','MarkerSize',25);
hp(2)=plotV(V(bcPrescribeList2,:),'g.','MarkerSize',25);
hp(3)=plotV(V(bcPrescribeList3,:),'b.','MarkerSize',25);
hp(4)=plotV(V(bcPrescribeList4,:),'y.','MarkerSize',25);
legend(hp,{'Node set 1','Node set 2','Node set 3','Node set 4'})
axisGeom(gca,fontSize);
camlight headlight;
drawnow; 

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='4.0'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

%Create control structure for use by all steps
stepStruct.Control.time_steps=numTimeSteps;
stepStruct.Control.step_size=1/numTimeSteps;
stepStruct.Control.solver.max_refs=max_refs;
stepStruct.Control.time_stepper.dtmin=dtmin;
stepStruct.Control.time_stepper.dtmax=dtmax; 
stepStruct.Control.time_stepper.max_retries=max_retries;
stepStruct.Control.time_stepper.opt_iter=opt_iter;

%Add template based default settings to proposed control section
[stepStruct.Control]=structComplete(stepStruct.Control,febio_spec.Control,1); %Complement provided with default if missing

%Remove control field (part of template) since step specific control sections are used
febio_spec=rmfield(febio_spec,'Control'); 

febio_spec.Step.step{1}.Control=stepStruct.Control;
febio_spec.Step.step{1}.ATTR.id=1;
febio_spec.Step.step{2}.Control=stepStruct.Control;
febio_spec.Step.step{2}.ATTR.id=2;

%Material section
materialName1='Material1';
febio_spec.Material.material{1}.ATTR.name=materialName1;
febio_spec.Material.material{1}.ATTR.type='Ogden';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c1=c1;
febio_spec.Material.material{1}.m1=m1;
febio_spec.Material.material{1}.c2=c1;
febio_spec.Material.material{1}.m2=-m1;
febio_spec.Material.material{1}.k=k;

% Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='penta6'; %Element type
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix
 
% -> NodeSets
nodeSetName1='bcPrescribeList1';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcPrescribeList1);

nodeSetName2='bcPrescribeList2';
febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.VAL=mrow(bcPrescribeList2);

nodeSetName3='bcPrescribeList3';
febio_spec.Mesh.NodeSet{3}.ATTR.name=nodeSetName3;
febio_spec.Mesh.NodeSet{3}.VAL=mrow(bcPrescribeList3);

nodeSetName4='bcPrescribeList4';
febio_spec.Mesh.NodeSet{4}.ATTR.name=nodeSetName4;
febio_spec.Mesh.NodeSet{4}.VAL=mrow(bcPrescribeList4);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

%Boundary condition section 

%STEP 1: Clamping compression -------------------------------------------
%Set 1

% Set 1 compression
febio_spec.Step.step{1}.Boundary.bc{1}.ATTR.name='bcPrescribeList01_01';
febio_spec.Step.step{1}.Boundary.bc{1}.ATTR.type='prescribed displacement';
febio_spec.Step.step{1}.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Step.step{1}.Boundary.bc{1}.dof='y';
febio_spec.Step.step{1}.Boundary.bc{1}.value.ATTR.lc=1;
febio_spec.Step.step{1}.Boundary.bc{1}.value.VAL=clampCompressiveDisplacement;
febio_spec.Step.step{1}.Boundary.bc{1}.relative=1;

febio_spec.Step.step{1}.Boundary.bc{2}.ATTR.name='FixedDisplacement01_01';
febio_spec.Step.step{1}.Boundary.bc{2}.ATTR.type='zero displacement';
febio_spec.Step.step{1}.Boundary.bc{2}.ATTR.node_set=nodeSetName1;
febio_spec.Step.step{1}.Boundary.bc{2}.x_dof=1;
febio_spec.Step.step{1}.Boundary.bc{2}.y_dof=0;
febio_spec.Step.step{1}.Boundary.bc{2}.z_dof=1;

% Set 2 compression
febio_spec.Step.step{1}.Boundary.bc{3}.ATTR.name='bcPrescribeList02_01';
febio_spec.Step.step{1}.Boundary.bc{3}.ATTR.type='prescribed displacement';
febio_spec.Step.step{1}.Boundary.bc{3}.ATTR.node_set=nodeSetName2;
febio_spec.Step.step{1}.Boundary.bc{3}.dof='y';
febio_spec.Step.step{1}.Boundary.bc{3}.value.ATTR.lc=1;
febio_spec.Step.step{1}.Boundary.bc{3}.value.VAL=-clampCompressiveDisplacement;
febio_spec.Step.step{1}.Boundary.bc{3}.relative=1;

febio_spec.Step.step{1}.Boundary.bc{4}.ATTR.name='FixedDisplacement02_01';
febio_spec.Step.step{1}.Boundary.bc{4}.ATTR.type='zero displacement';
febio_spec.Step.step{1}.Boundary.bc{4}.ATTR.node_set=nodeSetName2;
febio_spec.Step.step{1}.Boundary.bc{4}.x_dof=1;
febio_spec.Step.step{1}.Boundary.bc{4}.y_dof=0;
febio_spec.Step.step{1}.Boundary.bc{4}.z_dof=1;

% Set 3 compression
febio_spec.Step.step{1}.Boundary.bc{5}.ATTR.name='bcPrescribeList03_01';
febio_spec.Step.step{1}.Boundary.bc{5}.ATTR.type='prescribed displacement';
febio_spec.Step.step{1}.Boundary.bc{5}.ATTR.node_set=nodeSetName3;
febio_spec.Step.step{1}.Boundary.bc{5}.dof='y';
febio_spec.Step.step{1}.Boundary.bc{5}.value.ATTR.lc=1;
febio_spec.Step.step{1}.Boundary.bc{5}.value.VAL=-clampCompressiveDisplacement;
febio_spec.Step.step{1}.Boundary.bc{5}.relative=1;

febio_spec.Step.step{1}.Boundary.bc{6}.ATTR.name='FixedDisplacement03_01';
febio_spec.Step.step{1}.Boundary.bc{6}.ATTR.type='zero displacement';
febio_spec.Step.step{1}.Boundary.bc{6}.ATTR.node_set=nodeSetName3;
febio_spec.Step.step{1}.Boundary.bc{6}.x_dof=1;
febio_spec.Step.step{1}.Boundary.bc{6}.y_dof=0;
febio_spec.Step.step{1}.Boundary.bc{6}.z_dof=1;

% Set 4 compression
febio_spec.Step.step{1}.Boundary.bc{7}.ATTR.name='bcPrescribeList04_01';
febio_spec.Step.step{1}.Boundary.bc{7}.ATTR.type='prescribed displacement';
febio_spec.Step.step{1}.Boundary.bc{7}.ATTR.node_set=nodeSetName4;
febio_spec.Step.step{1}.Boundary.bc{7}.dof='y';
febio_spec.Step.step{1}.Boundary.bc{7}.value.ATTR.lc=1;
febio_spec.Step.step{1}.Boundary.bc{7}.value.VAL=clampCompressiveDisplacement;
febio_spec.Step.step{1}.Boundary.bc{7}.relative=1;

febio_spec.Step.step{1}.Boundary.bc{8}.ATTR.name='FixedDisplacement04_01';
febio_spec.Step.step{1}.Boundary.bc{8}.ATTR.type='zero displacement';
febio_spec.Step.step{1}.Boundary.bc{8}.ATTR.node_set=nodeSetName4;
febio_spec.Step.step{1}.Boundary.bc{8}.x_dof=1;
febio_spec.Step.step{1}.Boundary.bc{8}.y_dof=0;
febio_spec.Step.step{1}.Boundary.bc{8}.z_dof=1;

%STEP 2 Tension -------------------------------------------

% Set 1 tension
febio_spec.Step.step{2}.Boundary.bc{1}.ATTR.name='bcPrescribeList01_02';
febio_spec.Step.step{2}.Boundary.bc{1}.ATTR.type='prescribed displacement';
febio_spec.Step.step{2}.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Step.step{2}.Boundary.bc{1}.dof='x';
febio_spec.Step.step{2}.Boundary.bc{1}.value.ATTR.lc=1;
febio_spec.Step.step{2}.Boundary.bc{1}.value.VAL=0;
febio_spec.Step.step{2}.Boundary.bc{1}.relative=1;

febio_spec.Step.step{2}.Boundary.bc{2}.ATTR.name='bcPrescribeList02_02';
febio_spec.Step.step{2}.Boundary.bc{2}.ATTR.type='prescribed displacement';
febio_spec.Step.step{2}.Boundary.bc{2}.ATTR.node_set=nodeSetName1;
febio_spec.Step.step{2}.Boundary.bc{2}.dof='y';
febio_spec.Step.step{2}.Boundary.bc{2}.value.ATTR.lc=1;
febio_spec.Step.step{2}.Boundary.bc{2}.value.VAL=0;
febio_spec.Step.step{2}.Boundary.bc{2}.relative=1;

febio_spec.Step.step{2}.Boundary.bc{3}.ATTR.name='bcPrescribeList03_02';
febio_spec.Step.step{2}.Boundary.bc{3}.ATTR.type='prescribed displacement';
febio_spec.Step.step{2}.Boundary.bc{3}.ATTR.node_set=nodeSetName1;
febio_spec.Step.step{2}.Boundary.bc{3}.dof='z';
febio_spec.Step.step{2}.Boundary.bc{3}.value.ATTR.lc=2;
febio_spec.Step.step{2}.Boundary.bc{3}.value.VAL=clampTensionDisplacement;
febio_spec.Step.step{2}.Boundary.bc{3}.relative=1;

% Set 2 tension
febio_spec.Step.step{2}.Boundary.bc{4}.ATTR.name='bcPrescribeList04_02';
febio_spec.Step.step{2}.Boundary.bc{4}.ATTR.type='prescribed displacement';
febio_spec.Step.step{2}.Boundary.bc{4}.ATTR.node_set=nodeSetName2;
febio_spec.Step.step{2}.Boundary.bc{4}.dof='x';
febio_spec.Step.step{2}.Boundary.bc{4}.value.ATTR.lc=1;
febio_spec.Step.step{2}.Boundary.bc{4}.value.VAL=0;
febio_spec.Step.step{2}.Boundary.bc{4}.relative=1;

febio_spec.Step.step{2}.Boundary.bc{5}.ATTR.name='bcPrescribeList05_02';
febio_spec.Step.step{2}.Boundary.bc{5}.ATTR.type='prescribed displacement';
febio_spec.Step.step{2}.Boundary.bc{5}.ATTR.node_set=nodeSetName2;
febio_spec.Step.step{2}.Boundary.bc{5}.dof='y';
febio_spec.Step.step{2}.Boundary.bc{5}.value.ATTR.lc=1;
febio_spec.Step.step{2}.Boundary.bc{5}.value.VAL=0;
febio_spec.Step.step{2}.Boundary.bc{5}.relative=1;

febio_spec.Step.step{2}.Boundary.bc{6}.ATTR.name='bcPrescribeList06_02';
febio_spec.Step.step{2}.Boundary.bc{6}.ATTR.type='prescribed displacement';
febio_spec.Step.step{2}.Boundary.bc{6}.ATTR.node_set=nodeSetName2;
febio_spec.Step.step{2}.Boundary.bc{6}.dof='z';
febio_spec.Step.step{2}.Boundary.bc{6}.value.ATTR.lc=2;
febio_spec.Step.step{2}.Boundary.bc{6}.value.VAL=clampTensionDisplacement;
febio_spec.Step.step{2}.Boundary.bc{6}.relative=1;

% Set 3 tension
febio_spec.Step.step{2}.Boundary.bc{7}.ATTR.name='bcPrescribeList07_02';
febio_spec.Step.step{2}.Boundary.bc{7}.ATTR.type='prescribed displacement';
febio_spec.Step.step{2}.Boundary.bc{7}.ATTR.node_set=nodeSetName3;
febio_spec.Step.step{2}.Boundary.bc{7}.dof='x';
febio_spec.Step.step{2}.Boundary.bc{7}.value.ATTR.lc=1;
febio_spec.Step.step{2}.Boundary.bc{7}.value.VAL=0;
febio_spec.Step.step{2}.Boundary.bc{7}.relative=1;

febio_spec.Step.step{2}.Boundary.bc{8}.ATTR.name='bcPrescribeList08_02';
febio_spec.Step.step{2}.Boundary.bc{8}.ATTR.type='prescribed displacement';
febio_spec.Step.step{2}.Boundary.bc{8}.ATTR.node_set=nodeSetName3;
febio_spec.Step.step{2}.Boundary.bc{8}.dof='y';
febio_spec.Step.step{2}.Boundary.bc{8}.value.ATTR.lc=1;
febio_spec.Step.step{2}.Boundary.bc{8}.value.VAL=0;
febio_spec.Step.step{2}.Boundary.bc{8}.relative=1;

febio_spec.Step.step{2}.Boundary.bc{9}.ATTR.name='bcPrescribeList09_02';
febio_spec.Step.step{2}.Boundary.bc{9}.ATTR.type='prescribed displacement';
febio_spec.Step.step{2}.Boundary.bc{9}.ATTR.node_set=nodeSetName3;
febio_spec.Step.step{2}.Boundary.bc{9}.dof='z';
febio_spec.Step.step{2}.Boundary.bc{9}.value.ATTR.lc=1;
febio_spec.Step.step{2}.Boundary.bc{9}.value.VAL=0;
febio_spec.Step.step{2}.Boundary.bc{9}.relative=1;

% Set 4 tension
febio_spec.Step.step{2}.Boundary.bc{10}.ATTR.name='bcPrescribeList10_02';
febio_spec.Step.step{2}.Boundary.bc{10}.ATTR.type='prescribed displacement';
febio_spec.Step.step{2}.Boundary.bc{10}.ATTR.node_set=nodeSetName4;
febio_spec.Step.step{2}.Boundary.bc{10}.dof='x';
febio_spec.Step.step{2}.Boundary.bc{10}.value.ATTR.lc=1;
febio_spec.Step.step{2}.Boundary.bc{10}.value.VAL=0;
febio_spec.Step.step{2}.Boundary.bc{10}.relative=1;

febio_spec.Step.step{2}.Boundary.bc{11}.ATTR.name='bcPrescribeList11_02';
febio_spec.Step.step{2}.Boundary.bc{11}.ATTR.type='prescribed displacement';
febio_spec.Step.step{2}.Boundary.bc{11}.ATTR.node_set=nodeSetName4;
febio_spec.Step.step{2}.Boundary.bc{11}.dof='y';
febio_spec.Step.step{2}.Boundary.bc{11}.value.ATTR.lc=1;
febio_spec.Step.step{2}.Boundary.bc{11}.value.VAL=0;
febio_spec.Step.step{2}.Boundary.bc{11}.relative=1;

febio_spec.Step.step{2}.Boundary.bc{12}.ATTR.name='bcPrescribeList12_02';
febio_spec.Step.step{2}.Boundary.bc{12}.ATTR.type='prescribed displacement';
febio_spec.Step.step{2}.Boundary.bc{12}.ATTR.node_set=nodeSetName4;
febio_spec.Step.step{2}.Boundary.bc{12}.dof='z';
febio_spec.Step.step{2}.Boundary.bc{12}.value.ATTR.lc=1;
febio_spec.Step.step{2}.Boundary.bc{12}.value.VAL=0;
febio_spec.Step.step{2}.Boundary.bc{12}.relative=1;

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.name='LC_1';
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
%febio_spec.LoadData.load_controller{1}.extend='CONSTANT';
febio_spec.LoadData.load_controller{1}.points.pt.VAL=[0 0; 1 1; 2 1];

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{2}.ATTR.name='LC_2';
febio_spec.LoadData.load_controller{2}.ATTR.id=2;
febio_spec.LoadData.load_controller{2}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate='LINEAR';
%febio_spec.LoadData.load_controller{2}.extend='CONSTANT';
febio_spec.LoadData.load_controller{2}.points.pt.VAL=[0 0; 1 0; 2 1];

%Output section
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_strain;
febio_spec.Output.logfile.element_data{1}.ATTR.data='Ez';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

% Plotfile section
febio_spec.Output.plotfile.compression=0;

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
febioAnalysis.runMode=runMode;
febioAnalysis.maxLogCheckTime=120; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results 

if runFlag==1 %i.e. a succesful run
    
    %%     
    % Importing nodal displacements from a log file

    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),0,1);
    
    %Access data
    N_disp_mat=dataStruct.data; %Displacement
    timeVec=dataStruct.time; %Time
    
    %Create deformed coordinate set
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);

    
    %%
    % Importing element strain data from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_strain),0,1);
    
    %Access data
    E_strain=dataStruct.data;
    
    [F,CF]=element2patch(E,E_strain(:,:,end));
    CV=faceToVertexMeasure(F,V,CF);

    E_strain_middle_mean=squeeze(mean(E_strain(logicMiddleElements,:,:),1));

    %% Compute grimp implied strain
    gripStrainLinear=double(timeVec>1).*(timeVec-1).*appliedLinearStrain; 
    gripStrain=1/2*((gripStrainLinear+1).^2-1); %Green-Lagrange strain
    maxGridStrain=max(abs(gripStrain));

      %%
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations
        
    indBc=[bcPrescribeList1;bcPrescribeList2;bcPrescribeList3;bcPrescribeList4;];
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; hold on;%Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('$E_{zz}$','Interpreter','Latex')
    hp1=gpatch(Fb,V_DEF(:,:,end),CV,'k',1);
    hp1(1).FaceColor='interp';
    hp1(2).FaceColor='interp';
    hp2=plotV(V_DEF(indBc,:,end),'k.','MarkerSize',25);
    
    axisGeom(gca,fontSize);
    colormap(warmcold(250)); hc=colorbar;
    caxis([-maxGridStrain maxGridStrain]);
    hc.Ticks=linspace(-maxGridStrain,maxGridStrain,7);
    axis(axisLim(V_DEF)); %Set axis limits statically    

    camlight headlight; axis off;

    % Set up animation features
    animStruct.Time=timeVec; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments                
        CV=faceToVertexMeasure(E,V,E_strain(:,:,qt));

        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1(1) hp1(1) hp1(2) hp1(2) hp2 hp2 hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices','CData','XData','YData','ZData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),CV,V_DEF(:,:,qt),CV,V_DEF(indBc,1,qt),V_DEF(indBc,2,qt),V_DEF(indBc,3,qt)}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct); %Initiate animation feature
    gdrawnow;

    %%
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations
    
    indBc=[bcPrescribeList1;bcPrescribeList2;bcPrescribeList3;bcPrescribeList4;];
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; 
    subplot(1,2,1); hold on;
    title('$E_{zz}$','Interpreter','Latex');
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1=gpatch(Fb,V_DEF(:,:,end),CV,'none',1,0.5); 
    hp1(1).FaceColor='interp';
    hp1(2).FaceColor='interp';
    hp2=plotV(V_DEF(indBc,:,end),'k.','MarkerSize',1);    
    
    colormap(warmcold(250)); hc=colorbar;
    caxis([-maxGridStrain maxGridStrain]);
    hc.Ticks=linspace(-maxGridStrain,maxGridStrain,7);    
    axisGeom(gca,fontSize); camlight headlight; 
    axis(axisLim(V_DEF)); %Set axis limits statically
    view(0,0); axis off;
   
    subplot(1,2,2); hold on;    
    xlabel('Time (s)'); ylabel('$E_{zz}$','Interpreter','Latex');
    hpl4=plot(timeVec,E_strain_middle_mean,'g.-','LineWidth',3);
    hpl5=plot(timeVec,gripStrain,'r.-','LineWidth',3);
    hp4=plot(timeVec(end),E_strain_middle_mean(end),'g.','MarkerSize',50);
    hp5=plot(timeVec(end),gripStrain(end),'r.','MarkerSize',50);
    legend([hpl4 hpl5],{'True strain $E_{zz}$','"Intended applied" strain $E_{zz}$'},'Location','NorthOutside','Interpreter','Latex');
    axis tight; box on; grid on; set(gca,'FontSize',fontSize);

    % Set up animation features
    animStruct.Time=timeVec; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments
        V_def=V_DEF(:,:,qt); %Current nodal coordinates
        
        CV=faceToVertexMeasure(E,V,E_strain(:,:,qt));
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1(1) hp1(1) hp1(2) hp1(2) hp2 hp2 hp2 hp4 hp4 hp5 hp5]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices','CData','XData','YData','ZData','XData','YData','XData','YData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,CV,V_def,CV,V_def(indBc,1),V_def(indBc,2),V_def(indBc,3),timeVec(qt),E_strain_middle_mean(qt),timeVec(qt),gripStrain(qt)}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct); %Initiate animation feature
    drawnow;
       
end

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
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
