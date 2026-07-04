%% DEMO_febio_0021_bone_three_point_bending
% Below is a demonstration for:
%
% * Building geometry for a bone
% * Defining the boundary conditions
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement results

%% Keywords
%
% * febio_spec version 4.0
% * febio, FEBio
% * indentation
% * contact, sliding, sticky, friction
% * rigid body constraints
% * tetrahedral elements, tet4
% * triangular elements, tri3
% * three point bending
% * static, solid
% * hyperelastic, Ogden
% * displacement logfile
% * stress logfile

%%

clear; close all; clc;

%% Plot settings
fontSize=15;
faceAlpha1=0.8;
faceAlpha2=0.3;
markerSize=40;
lineWidth=3;

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');
pathNameSTL=fullfile(defaultFolder,'data','STL'); 

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress

%Geometric parameters
corticalThickness=3; %Thickness used for cortical material definition
volumeFactor=10; %Factor to scale desired volume for interior elements w.r.t. boundary elements
pointSpacing=5;

%Material parameter set
D_density=1e-9; %Density for all materials

E_youngs1=17000; %Youngs modulus
nu1=0.25; %Poissons ratio

% Cancellous bone
E_youngs2=1500; %Youngs modulus
nu2=0.25; %Poissons ratio

% FEA control settings
analysisType='STATIC';
numTimeSteps=20; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=15; %Optimum number of iterations
max_retries=6; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
min_residual=1e-10;
symmetric_stiffness=0;

runMode='external';% 'internal' or 'external'

%Contact parameters
contactInitialOffset=0.1;
contactPenalty=1000;
laugon=0;
minaug=1;
maxaug=10;

%Boundary conditions specification
bcFix=0;

%Loading parameters
zDisp=-20-2*contactInitialOffset; 

%% Get bone geometry

%Load basic model
[F_bone,V_bone]=graphicsModels(5); 
V_bone=V_bone*1000; %Scale to mm

%Reorient
V_mean=mean(V_bone,1);
V_bone=V_bone-V_mean(ones(size(V_bone,1),1),:); %Center around origin
[R]=pointSetPrincipalDir(V_bone); %Get rotation matrix
V_bone=V_bone*R; %Rotate

%% Remesh using geomgram

optionStruct.pointSpacing=pointSpacing;
optionStruct.disp_on=0;
[F_bone,V_bone]=ggremesh(F_bone,V_bone,optionStruct);

%% Create and position cylinder geometry

pointSpacingBeams=pointSpacing/2; 

inputStruct.cylRadius=10;
inputStruct.numRadial=round((2*pi*inputStruct.cylRadius)./pointSpacingBeams);
inputStruct.cylHeight=max(V_bone(:,2))-min(V_bone(:,2));
nh=round(inputStruct.cylHeight./pointSpacingBeams);
nh=nh+double(iseven(nh));
inputStruct.numHeight=nh;
inputStruct.meshType='quad';
inputStruct.closeOpt=0;

% Derive patch data for a cylinder
[Fc,Vc,Cc]=patchcylinder(inputStruct);
R=euler2DCM([0.5*pi 0 0]);
Vc=Vc*R;
Vc1=Vc;
Vc2=Vc;
Vc2(:,1)=Vc2(:,1)-60;
Vc3=Vc;
Vc3(:,1)=Vc3(:,1)+60;

logicSelect=min(Vc1(:,1))<V_bone(:,1) & max(Vc1(:,1))>V_bone(:,1);
zOffset=max(V_bone(logicSelect,3));
Vc1(:,3)=Vc1(:,3)-min(Vc1(:,3))+zOffset+contactInitialOffset;

logicSelect=min(Vc2(:,1))<V_bone(:,1) & max(Vc2(:,1))>V_bone(:,1);
zOffset=min(V_bone(logicSelect,3));
Vc2(:,3)=Vc2(:,3)-max(Vc2(:,3))+zOffset-contactInitialOffset;

logicSelect=min(Vc3(:,1))<V_bone(:,1) & max(Vc3(:,1))>V_bone(:,1);
zOffset=min(V_bone(logicSelect,3));
Vc3(:,3)=Vc3(:,3)-max(Vc3(:,3))+zOffset-contactInitialOffset;

%%
% Plotting surface geometry

cFigure;
hold on;
gpatch(F_bone,V_bone,'kw','k',1);
gpatch(Fc,Vc1,'gw','g',1);
gpatch(Fc,Vc2,'rw','r',1);
gpatch(Fc,Vc3,'bw','b',1);

axisGeom;
camlight headlight;
drawnow

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
tetGenStruct.faceBoundaryMarker=ones(size(F_bone,1),1); %Face boundary markers
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

%% Joining node sets

Fc1=Fc+size(V,1);
Fc2=Fc+size(V,1)+size(Vc1,1);
Fc3=Fc+size(V,1)+size(Vc1,1)+size(Vc2,1);
V=[V;Vc1;Vc2;Vc3];

%% Define boundary conditions

if bcFix==1
    %Supported nodes
    indb=unique(Fb(:));
    logicLeft=V(indb,1)<-200;
    bcSupportList=indb(logicLeft);
    
    %%
    % Visualize BC's
    hf=cFigure;
    title('Boundary conditions model','FontSize',fontSize);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
    hold on;
    
    gpatch(Fb,V,'kw','none',faceAlpha2);
    hp1(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',markerSize);
    
    legend(hp1,{'BC support'});
    
    axisGeom(gca,fontSize);
    camlight headlight;
    drawnow;
    
end

%% Define bone contact surfaces

N=patchNormal(fliplr(Fb),V);
nz=[0 0 -1];
d=dot(N,nz(ones(size(N,1),1),:),2);

contactAdd=pointSpacing;

VF=patchCentre(Fb,V);
logicSecondary1=VF(:,1)>(min(Vc1(:,1))-contactAdd) & VF(:,1)<(max(Vc1(:,1))+contactAdd) & d<0;
logicSecondary2=VF(:,1)>(min(Vc2(:,1))-contactAdd) & VF(:,1)<(max(Vc2(:,1))+contactAdd) & d>0;
logicSecondary3=VF(:,1)>(min(Vc3(:,1))-contactAdd) & VF(:,1)<(max(Vc3(:,1))+contactAdd) & d>0;

F_prim1=fliplr(Fb(logicSecondary1,:));
F_prim2=fliplr(Fb(logicSecondary2,:));
F_prim3=fliplr(Fb(logicSecondary3,:));

%%
% Visualize

cFigure; hold on;

gpatch(Fb,V,'w','none',0.5);

hp2(1)=gpatch(F_prim1,V,'r','k',1);
patchNormPlot(F_prim1,V);

hp2(2)=gpatch(F_prim2,V,'g','k',1);
patchNormPlot(F_prim2,V);

hp2(3)=gpatch(F_prim3,V,'b','k',1);
patchNormPlot(F_prim3,V);

hp2(4)=gpatch(Fc1,V,'c','k',1);
patchNormPlot(Fc1,V);

hp2(5)=gpatch(Fc2,V,'y','k',1);
patchNormPlot(Fc2,V);

hp2(6)=gpatch(Fc3,V,'m','k',1);
patchNormPlot(Fc3,V);

legend(hp2,{'primary 1','primary 2','primary 3','secondary 1','secondary 2','secondary 3'});

axisGeom;
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

%Control section
febio_spec.Control.analysis='STATIC';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.solver.qn_method.max_ups=max_ups;
febio_spec.Control.solver.symmetric_stiffness=symmetric_stiffness;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;

%Material section
materialName1='Material1';
febio_spec.Material.material{1}.ATTR.name=materialName1;
febio_spec.Material.material{1}.ATTR.type='neo-Hookean';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.E=E_youngs1;
febio_spec.Material.material{1}.v=nu1;
febio_spec.Material.material{1}.density=D_density;

materialName2='Material2';
febio_spec.Material.material{2}.ATTR.name=materialName2;
febio_spec.Material.material{2}.ATTR.type='neo-Hookean';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.E=E_youngs2;
febio_spec.Material.material{2}.v=nu2;
febio_spec.Material.material{2}.density=D_density;

materialName3='Material3';
febio_spec.Material.material{3}.ATTR.name=materialName3;
febio_spec.Material.material{3}.ATTR.type='rigid body';
febio_spec.Material.material{3}.ATTR.id=3;
febio_spec.Material.material{3}.density=D_density;
febio_spec.Material.material{3}.center_of_mass=mean(Vc1,1);

materialName4='Material4';
febio_spec.Material.material{4}.ATTR.name=materialName4;
febio_spec.Material.material{4}.ATTR.type='rigid body';
febio_spec.Material.material{4}.ATTR.id=4;
febio_spec.Material.material{4}.density=D_density;
febio_spec.Material.material{4}.center_of_mass=mean(Vc2,1);

materialName5='Material5';
febio_spec.Material.material{5}.ATTR.name=materialName5;
febio_spec.Material.material{5}.ATTR.type='rigid body';
febio_spec.Material.material{5}.ATTR.id=5;
febio_spec.Material.material{5}.density=D_density;
febio_spec.Material.material{5}.center_of_mass=mean(Vc3,1);

%Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='CorticalBone';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='tet4'; %Element type of this set
febio_spec.Mesh.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E1,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E1;

partName2='CancellousBone';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='tet4'; %Element type of this set
febio_spec.Mesh.Elements{2}.ATTR.mat=2; %material index for this set 
febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=E2;

partName3='Part3';
febio_spec.Mesh.Elements{3}.ATTR.name=partName3; %Name of this part
febio_spec.Mesh.Elements{3}.ATTR.type='quad4'; %Element type of this set
febio_spec.Mesh.Elements{3}.ATTR.mat=3; %material index for this set
febio_spec.Mesh.Elements{3}.elem.ATTR.id=size(E1,1)+size(E2,1)+(1:1:size(Fc1,1))'; %Element id's
febio_spec.Mesh.Elements{3}.elem.VAL=Fc1;

partName4='Part4';
febio_spec.Mesh.Elements{4}.ATTR.name=partName4; %Name of this part
febio_spec.Mesh.Elements{4}.ATTR.type='quad4'; %Element type of this set
febio_spec.Mesh.Elements{4}.ATTR.mat=4; %material index for this set
febio_spec.Mesh.Elements{4}.elem.ATTR.id=size(E1,1)+size(E2,1)+size(Fc1,1)+(1:1:size(Fc2,1))'; %Element id's
febio_spec.Mesh.Elements{4}.elem.VAL=Fc2;

partName5='Part5';
febio_spec.Mesh.Elements{5}.ATTR.name=partName5; %Name of this part
febio_spec.Mesh.Elements{5}.ATTR.type='quad4'; %Element type of this set
febio_spec.Mesh.Elements{5}.ATTR.mat=5; %material index for this set
febio_spec.Mesh.Elements{5}.elem.ATTR.id=size(E1,1)+size(E2,1)+size(Fc1,1)+size(Fc2,1)+(1:1:size(Fc3,1))'; %Element id's
febio_spec.Mesh.Elements{5}.elem.VAL=Fc3;

%MeshDomains section
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat=materialName1;

febio_spec.MeshDomains.SolidDomain{2}.ATTR.name=partName2;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat=materialName2;

febio_spec.MeshDomains.ShellDomain{1}.ATTR.name=partName3;
febio_spec.MeshDomains.ShellDomain{1}.ATTR.mat=materialName3;

febio_spec.MeshDomains.ShellDomain{2}.ATTR.name=partName4;
febio_spec.MeshDomains.ShellDomain{2}.ATTR.mat=materialName4;

febio_spec.MeshDomains.ShellDomain{3}.ATTR.name=partName5;
febio_spec.MeshDomains.ShellDomain{3}.ATTR.mat=materialName5;

% -> NodeSets
if bcFix==1
    nodeSetName1='bcSupportList';
    febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
    febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcSupportList);
end

% -> Surfaces
surfaceNameSec1='contact_secondary1';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceNameSec1;
febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(Fc1,1))';
febio_spec.Mesh.Surface{1}.quad4.VAL=Fc1;

surfaceNameSec2='contact_secondary2';
febio_spec.Mesh.Surface{2}.ATTR.name=surfaceNameSec2;
febio_spec.Mesh.Surface{2}.quad4.ATTR.id=(1:1:size(Fc2,1))';
febio_spec.Mesh.Surface{2}.quad4.VAL=Fc2;

surfaceNameSec3='contact_secondary3';
febio_spec.Mesh.Surface{3}.ATTR.name=surfaceNameSec3;
febio_spec.Mesh.Surface{3}.quad4.ATTR.id=(1:1:size(Fc3,1))';
febio_spec.Mesh.Surface{3}.quad4.VAL=Fc3;

surfaceNamePrim1='contact_primary1';
febio_spec.Mesh.Surface{4}.ATTR.name=surfaceNamePrim1;
febio_spec.Mesh.Surface{4}.tri3.ATTR.id=(1:1:size(F_prim1,1))';
febio_spec.Mesh.Surface{4}.tri3.VAL=F_prim1;

surfaceNamePrim2='contact_primary2';
febio_spec.Mesh.Surface{5}.ATTR.name=surfaceNamePrim2;
febio_spec.Mesh.Surface{5}.tri3.ATTR.id=(1:1:size(F_prim2,1))';
febio_spec.Mesh.Surface{5}.tri3.VAL=F_prim2;

surfaceNamePrim3='contact_primary3';
febio_spec.Mesh.Surface{6}.ATTR.name=surfaceNamePrim3;
febio_spec.Mesh.Surface{6}.tri3.ATTR.id=(1:1:size(F_prim3,1))';
febio_spec.Mesh.Surface{6}.tri3.VAL=F_prim3;

% -> Surface pairs
febio_spec.Mesh.SurfacePair{1}.ATTR.name='Contact1';
febio_spec.Mesh.SurfacePair{1}.primary=surfaceNamePrim1;
febio_spec.Mesh.SurfacePair{1}.secondary=surfaceNameSec1;

febio_spec.Mesh.SurfacePair{2}.ATTR.name='Contact2';
febio_spec.Mesh.SurfacePair{2}.primary=surfaceNamePrim2;
febio_spec.Mesh.SurfacePair{2}.secondary=surfaceNameSec2;

febio_spec.Mesh.SurfacePair{3}.ATTR.name='Contact3';
febio_spec.Mesh.SurfacePair{3}.primary=surfaceNamePrim3;
febio_spec.Mesh.SurfacePair{3}.secondary=surfaceNameSec3;

%Boundary condition section 
% -> Fix boundary conditions
if bcFix==1
    febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_xyz';
    febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
    febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
    febio_spec.Boundary.bc{1}.x_dof=1;
    febio_spec.Boundary.bc{1}.y_dof=1;
    febio_spec.Boundary.bc{1}.z_dof=1;
end

%Rigid section 
% ->Rigid body fix boundary conditions
febio_spec.Rigid.rigid_bc{1}.ATTR.name='RotFix_1';
febio_spec.Rigid.rigid_bc{1}.ATTR.type='rigid_fixed';
febio_spec.Rigid.rigid_bc{1}.rb=3;
febio_spec.Rigid.rigid_bc{1}.Rx_dof=1;
febio_spec.Rigid.rigid_bc{1}.Ry_dof=1;
%febio_spec.Rigid.rigid_bc{1}.Rz_dof=1;
febio_spec.Rigid.rigid_bc{1}.Ru_dof=1;
febio_spec.Rigid.rigid_bc{1}.Rv_dof=1;
febio_spec.Rigid.rigid_bc{1}.Rw_dof=1;

febio_spec.Rigid.rigid_bc{2}.ATTR.name='RotFix_2';
febio_spec.Rigid.rigid_bc{2}.ATTR.type='rigid_fixed';
febio_spec.Rigid.rigid_bc{2}.rb=4;
febio_spec.Rigid.rigid_bc{2}.Rx_dof=1;
febio_spec.Rigid.rigid_bc{2}.Ry_dof=1;
febio_spec.Rigid.rigid_bc{2}.Rz_dof=1;
febio_spec.Rigid.rigid_bc{2}.Ru_dof=1;
febio_spec.Rigid.rigid_bc{2}.Rv_dof=1;
febio_spec.Rigid.rigid_bc{2}.Rw_dof=1;

febio_spec.Rigid.rigid_bc{3}.ATTR.name='RotFix_3';
febio_spec.Rigid.rigid_bc{3}.ATTR.type='rigid_fixed';
febio_spec.Rigid.rigid_bc{3}.rb=5;
febio_spec.Rigid.rigid_bc{3}.Rx_dof=1;
febio_spec.Rigid.rigid_bc{3}.Ry_dof=1;
febio_spec.Rigid.rigid_bc{3}.Rz_dof=1;
febio_spec.Rigid.rigid_bc{3}.Ru_dof=1;
febio_spec.Rigid.rigid_bc{3}.Rv_dof=1;
febio_spec.Rigid.rigid_bc{3}.Rw_dof=1;

% ->Rigid body prescribe boundary conditions
febio_spec.Rigid.rigid_bc{4}.ATTR.name='RigidPrescribe';
febio_spec.Rigid.rigid_bc{4}.ATTR.type='rigid_displacement';
febio_spec.Rigid.rigid_bc{4}.rb=3;
febio_spec.Rigid.rigid_bc{4}.dof='z';
febio_spec.Rigid.rigid_bc{4}.value.ATTR.lc=1;
febio_spec.Rigid.rigid_bc{4}.value.VAL=zDisp;
febio_spec.Rigid.rigid_bc{4}.relative=0;

%Contact section
for qc=1:1:3
    febio_spec.Contact.contact{qc}.ATTR.surface_pair=febio_spec.Mesh.SurfacePair{qc}.ATTR.name;
    febio_spec.Contact.contact{qc}.ATTR.type='sticky';
    febio_spec.Contact.contact{qc}.penalty=contactPenalty;
    febio_spec.Contact.contact{qc}.laugon=laugon;
    febio_spec.Contact.contact{qc}.tolerance=0.2;
    febio_spec.Contact.contact{qc}.minaug=minaug;
    febio_spec.Contact.contact{qc}.maxaug=maxaug;
    febio_spec.Contact.contact{qc}.snap_tol=0.01;
    febio_spec.Contact.contact{qc}.max_traction=0;
    febio_spec.Contact.contact{qc}.search_tolerance=0.01;
end

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.name='LC_1';
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
%febio_spec.LoadData.load_controller{1}.extend='CONSTANT';
febio_spec.LoadData.load_controller{1}.points.pt.VAL=[0 0; 1 1];

%Output section
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.rigid_body_data{1}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.rigid_body_data{1}.ATTR.data='Fx;Fy;Fz';
febio_spec.Output.logfile.rigid_body_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.rigid_body_data{1}.VAL=3;

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='s1';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

%% Quick viewing of the FEBio input file structure
% The |febView| function can be used to view the xml structure in a MATLAB
% figure window. 

%%
% |febView(febio_spec); %Viewing the febio file|

%% Exporting the FEBio input file
% Exporting the febio_spec structure to an FEBio input file is done using
% the |febioStruct2xml| function. 

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
%system(['gedit ',febioFebFileName,' &']);

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
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

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
    % Importing element stress from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),0,1);
    
    %Access data
    E_stress_mat=dataStruct.data;
    E_stress_mat(isnan(E_stress_mat))=0;
    
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    [CV]=faceToVertexMeasure(E,V,E_stress_mat(:,:,end));
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('$\sigma_{1}$ [MPa]','Interpreter','Latex')
    hp=gpatch(Fb,V_DEF(:,:,end),CV,'k',1); %Add graphics object to animate
    hp.FaceColor='interp';
    hp2=gpatch([Fc1;Fc2;Fc3],V,'w','k',0.5); 
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([min(E_stress_mat(:)) max(E_stress_mat(:))]/20);    
    axis(axisLim(V_DEF)); %Set axis limits statically    
    camlight headlight;        
        
    % Set up animation features
    animStruct.Time=timeVec; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        
        [CV]=faceToVertexMeasure(E,V,E_stress_mat(:,:,qt));
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),CV,V_DEF(:,:,qt)}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;
    
    %%
    % Importing rigidbody force data from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_force),0,1);
    
    %Access data
    Force_mat=dataStruct.data;
    
    Fz=squeeze(Force_mat(:,3,:));
       
    %%    
    % Visualize stress-stretch curve
    
    cFigure; hold on;    
    title('Force-displacement curve','FontSize',fontSize);
    xlabel('Displacement [mm]','FontSize',fontSize,'Interpreter','Latex'); 
    ylabel('$F_z$ [N]','FontSize',fontSize,'Interpreter','Latex'); 
    
    plot(timeVec(:).*zDisp,Fz(:),'r-','lineWidth',lineWidth);
    
    view(2); axis tight;  grid on; axis square; box on; 
    set(gca,'FontSize',fontSize);
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
