%% DEMO_febio_0044_mammography_01.m
% Below is a demonstration for:
%
% * Building geometry for a slab with hexahedral elements, and a
% triangulated sphere.
% * Defining the boundary conditions
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement results

%% Keywords
%
% * febio_spec version 2.5
% * febio, FEBio
% * indentation
% * contact, sliding, sticky, friction
% * rigid body constraints
% * tetrahedral elements, tet8
% * triangular elements, tri3
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

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

%Breast geometry
r=40; %Breast radius
r1=r/2.5;
r2=r/7;
rm=mean([r1 r2]);
w=(r1-r2)/20;
h=r2;
dx=r/2; %Gravity direction shape alteration factor
nRefine=3; %Number of refine steps for hemi-sphere
plateDisplacement=25;
volumeFactor=3;
loadAngle=(0/180)*pi;

%Material parameter set
c1_1=1e-3; %Shear-modulus-like parameter
m1_1=6; %Material parameter setting degree of non-linearity
k_factor=1e2; %Bulk modulus factor
k_1=c1_1*k_factor; %Bulk modulus

c1_2=5e-3; %Shear-modulus-like parameter
m1_2=2; %Material parameter setting degree of non-linearity
k_factor=1e2; %Bulk modulus factor
k_2=c1_2*k_factor; %Bulk modulus

% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=10; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
symmetric_stiffness=0;
min_residual=1e-20;

%Contact parameters
contactInitialOffset=0.1;
contactAlg=5;
contactPenalty=10;
laugon=0;
minaug=1;
maxaug=10;
fric_coeff=0.1;
switch contactAlg
    case 1
        contactType='sticky';
    case 2
        contactType='facet-to-facet sliding';
    case 3
        contactType='sliding_with_gaps';
    case 4
        contactType='sliding2';
    case 5
        contactType='sliding-elastic';
end

%% Create hemi-sphere
[F,V,C_hemiSphereLabel]=hemiSphereMesh(nRefine,r,1); %Construct hemi-shere mesh
pointSpacing=mean(patchEdgeLengths(F,V)); % Get point spacing from mesh

% %%
% % Visualize hemi-sphere
% cFigure; hold on;
% gpatch(F,V,C_hemiSphereLabel);
% axisGeom;
% camlight headlight;
% icolorbar;
% gdrawnow;

%% Change shape of hemi-sphere to create basic breast model

indExclude=unique(F(C_hemiSphereLabel==2,:));
logicExclude=false(size(V,1),1);
logicExclude(indExclude)=1;

dt=sqrt(sum(V(:,[1 2]).^2,2));
logicHigh1=dt<r1 & ~logicExclude;
logicHigh2=dt<r2 & ~logicExclude;

C_skin=double(logicHigh1);
C_skin(logicHigh2)=2;

t=linspace(0,2*pi,500);
x=rm*sin(t);
y=rm*cos(t);
vc=[x(:) y(:)];

[d]=minDist(V(C_skin==1,[1 2]),vc);
dtt=d.^3;
dtt=dtt-min(dtt);
dtt=dtt./max(dtt);
dtt=abs(dtt-1)*w;
V(C_skin==1,3)=V(C_skin==1,3)+dtt;

f=V(:,3);
f=f-min(f(:));
f=f./max(f(:));
V(:,1)=V(:,1)+dx.*f;

dtt=dt(C_skin==2).^3;
dtt=dtt-min(dtt);
dtt=dtt./max(dtt);
dtt=abs(dtt-1)*h;
V(C_skin==2,3)=V(C_skin==2,3)+dtt;

%% Rotate model
R=euler2DCM([pi -0.5*pi 0]);
V=V*R;

%%

[Fs,Vs]=geoSphere(2,r/6);
Vs(:,1)=Vs(:,1)-r/2;
Vs(:,2)=Vs(:,2)-r/4;
Vs(:,3)=Vs(:,3)-r/2;

%%

C=[C_hemiSphereLabel;(max(C_hemiSphereLabel(:))+1)*ones(size(Fs,1),1)];
F=[F;Fs+size(V,1)];
V=[V;Vs];

cFigure; hold on;
gpatch(F,V,C,'none',0.5);
axisGeom;
camlight headlight;
icolorbar;
gdrawnow;

%% Get interior points
[V_in1]=getInnerPoint({F,Fs},{V,Vs});
[V_in2]=getInnerPoint(Fs,Vs);

V_regions=[V_in1; V_in2];

cFigure; hold on;
gpatch(F,V,C,'none',0.5);
% gpatch(Fs,Vs,'kw','none',1);
plotV(V_regions,'k.','markerSize',50);
axisGeom;
camlight headlight;
gdrawnow;

%%

faceBoundaryMarker=C;

[regionA]=tetVolMeanEst(F,V); %Volume for regular tets
inputStruct.stringOpt='-pq1.2AaY';
inputStruct.Faces=fliplr(F);
inputStruct.Nodes=V;
inputStruct.holePoints=[];
inputStruct.faceBoundaryMarker=faceBoundaryMarker; %Face boundary markers
inputStruct.regionPoints=V_regions; %region points
inputStruct.regionA=regionA*ones(size(V_regions,1),1)*volumeFactor;
inputStruct.minRegionMarker=2; %Minimum region marker

% Mesh model using tetrahedral elements using tetGen
[meshOutput]=runTetGen(inputStruct); %Run tetGen

% Access model element and patch data
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
V=meshOutput.nodes;
CE=meshOutput.elementMaterialID;
E=meshOutput.elements;

%%

Q=euler2DCM([loadAngle 0 0]);
V=V*Q; 

%%

cFigure; hold on;
hp=gpatch(Fb,V,Cb,'none',0.5);
axisGeom;
camlight headlight;
icolorbar;
gdrawnow;

%% Visualizing mesh using |meshView|, see also |anim8|

meshView(meshOutput);

%% Split element sets

E1=E(meshOutput.elementMaterialID==-2,:);
E2=E(meshOutput.elementMaterialID==-3,:);

%% Building plate models

%Basic side curve of plate
Vt=[0 0 0; 0 0 15; 1.5*r 0 15; 1.5*r 0 0; ];

%Fillet side curve of plate
rFillet=6; %Fillet radius
np=25; %Number of points used to construct each fillet edge
closedLoopOption=0; %Use 1 if curve represents a closed loop but containes unique points
[Vc]=filletCurve(Vt,rFillet,np,closedLoopOption);

%Extrude to form plate
cPar.pointSpacing=pointSpacing/2;
cPar.depth=3*r;
cPar.patchType='quad';
cPar.dir=0;
cPar.closeLoopOpt=0;
[Fp1,Vp1]=polyExtrude(Vc,cPar);
Fp1=fliplr(Fp1);
Vp1(:,3)=Vp1(:,3)-max(Vp1(:,3))+min(V(:,3));
Vp1(:,1)=Vp1(:,1)-max(Vp1(:,1))-r/7;

%Copy to create second plate
Fp2=fliplr(Fp1);
Vp2=Vp1;
Vp2(:,3)=-Vp2(:,3);
Vp2(:,3)=Vp2(:,3)-min(Vp2(:,3))+max(V(:,3));

cFigure; hold on;
gpatch(Fb,V,'gw','k',1);
gpatch(Fp1,Vp1,'kw','bw',1);
gpatch(Fp2,Vp2,'kw','rw',1);
axisGeom;
camlight headlight;
gdrawnow;

%% Joining node sets

Fp1=Fp1+size(V,1); %Fixed element indices
Fp2=Fp2+size(V,1)+size(Vp1,1); %Fixed element indices
V=[V;Vp1;Vp2;]; %Combined node sets

V=V*Q';
loadDir=[0 0 1]*Q';
plateDisplacement_XYZ=plateDisplacement.*loadDir;

%%
% Plotting joined geometry
cFigure;
title('Joined node sets','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
gpatch(Fb,V,'gw','k',1);
gpatch(Fp1,V,'kw','bw',1);
gpatch(Fp2,V,'kw','rw',1);
axisGeom;
camlight headlight;
gdrawnow;

%% Define contact surfaces

% The rigid master surface of the sphere
F_contact_master1=Fp1;
F_contact_master2=Fp2;

% The deformable slave surface of the slab
logicContactSurf1=Cb==1;
F_contact_slave=fliplr(Fb(logicContactSurf1,:));

% Plotting surface models
cFigure; hold on;
title('Contact sets and normal directions','FontSize',fontSize);

gpatch(Fb,V,'kw','none',faceAlpha2);
hl(1)=gpatch(F_contact_master1,V,'gw','k',1);
patchNormPlot(F_contact_master1,V);
hl(2)=gpatch(F_contact_master2,V,'rw','k',1);
patchNormPlot(F_contact_master2,V);
hl(3)=gpatch(F_contact_slave,V,'bw','k',1);
patchNormPlot(F_contact_slave,V);

legend(hl,{'Master 1','Master 2','Slave'});

axisGeom(gca,fontSize);
camlight headlight;
gdrawnow;

%% Define boundary conditions

%Supported nodes
logicRigid=Cb==2;
Fr=Fb(logicRigid,:);
bcSupportList=unique(Fr(:));

%%
% Visualize BC's
hf=cFigure;
title('Boundary conditions model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,V,'kw','none',faceAlpha2);

hl2(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',markerSize);

legend(hl2,{'BC support'});

axisGeom(gca,fontSize);
camlight headlight;
gdrawnow;

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
febio_spec.Control.symmetric_stiffness=symmetric_stiffness;
febio_spec.Control.min_residual=min_residual;

%Material section
febio_spec.Material.material{1}.ATTR.type='Ogden';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c1=c1_1;
febio_spec.Material.material{1}.m1=m1_1;
febio_spec.Material.material{1}.c2=c1_1;
febio_spec.Material.material{1}.m2=-m1_1;
febio_spec.Material.material{1}.k=k_1;

febio_spec.Material.material{2}.ATTR.type='Ogden';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.c1=c1_2;
febio_spec.Material.material{2}.m1=m1_2;
febio_spec.Material.material{2}.c2=c1_2;
febio_spec.Material.material{2}.m2=-m1_2;
febio_spec.Material.material{2}.k=k_2;

febio_spec.Material.material{3}.ATTR.type='rigid body';
febio_spec.Material.material{3}.ATTR.id=3;
febio_spec.Material.material{3}.density=1;
febio_spec.Material.material{3}.center_of_mass=mean(Vp1,1);

febio_spec.Material.material{4}.ATTR.type='rigid body';
febio_spec.Material.material{4}.ATTR.id=4;
febio_spec.Material.material{4}.density=1;
febio_spec.Material.material{4}.center_of_mass=mean(Vp2,1);

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='tet4'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set
febio_spec.Geometry.Elements{1}.ATTR.name='breastNormal'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(E1,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=E1;

febio_spec.Geometry.Elements{2}.ATTR.type='tet4'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat=2; %material index for this set
febio_spec.Geometry.Elements{2}.ATTR.name='breastInclusion'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL=E2;

febio_spec.Geometry.Elements{3}.ATTR.type='quad4'; %Element type of this set
febio_spec.Geometry.Elements{3}.ATTR.mat=3; %material index for this set
febio_spec.Geometry.Elements{3}.ATTR.name='PlateBottom'; %Name of the element set
febio_spec.Geometry.Elements{3}.elem.ATTR.id=size(E1,1)+size(E2,1)+(1:1:size(Fp1,1))'; %Element id's
febio_spec.Geometry.Elements{3}.elem.VAL=Fp1;

febio_spec.Geometry.Elements{4}.ATTR.type='quad4'; %Element type of this set
febio_spec.Geometry.Elements{4}.ATTR.mat=4; %material index for this set
febio_spec.Geometry.Elements{4}.ATTR.name='PlateTop'; %Name of the element set
febio_spec.Geometry.Elements{4}.elem.ATTR.id=size(E1,1)+size(E2,1)+size(Fp1,1)+(1:1:size(Fp2,1))'; %Element id's
febio_spec.Geometry.Elements{4}.elem.VAL=Fp2;

% -> NodeSets
febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcSupportList(:);

% -> Surfaces
febio_spec.Geometry.Surface{1}.ATTR.name='contact_master1';
febio_spec.Geometry.Surface{1}.quad4.ATTR.lid=(1:1:size(F_contact_master1,1))';
febio_spec.Geometry.Surface{1}.quad4.VAL=F_contact_master1;

febio_spec.Geometry.Surface{2}.ATTR.name='contact_master2';
febio_spec.Geometry.Surface{2}.quad4.ATTR.lid=(1:1:size(F_contact_master2,1))';
febio_spec.Geometry.Surface{2}.quad4.VAL=F_contact_master2;

febio_spec.Geometry.Surface{3}.ATTR.name='contact_slave';
febio_spec.Geometry.Surface{3}.tri3.ATTR.lid=(1:1:size(F_contact_slave,1))';
febio_spec.Geometry.Surface{3}.tri3.VAL=F_contact_slave;

% -> Surface pairs
febio_spec.Geometry.SurfacePair{1}.ATTR.name='Contact1';
febio_spec.Geometry.SurfacePair{1}.master.ATTR.surface=febio_spec.Geometry.Surface{1}.ATTR.name;
febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface=febio_spec.Geometry.Surface{3}.ATTR.name;

febio_spec.Geometry.SurfacePair{2}.ATTR.name='Contact2';
febio_spec.Geometry.SurfacePair{2}.master.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;
febio_spec.Geometry.SurfacePair{2}.slave.ATTR.surface=febio_spec.Geometry.Surface{3}.ATTR.name;

%Boundary condition section
% -> Fix boundary conditions
febio_spec.Boundary.fix{1}.ATTR.bc='x';
febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{2}.ATTR.bc='y';
febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{3}.ATTR.bc='z';
febio_spec.Boundary.fix{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;

% -> Prescribed boundary conditions on the rigid body
febio_spec.Boundary.rigid_body{1}.ATTR.mat=3;
febio_spec.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='Rx';
febio_spec.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='Ry';
febio_spec.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='Rz';

febio_spec.Boundary.rigid_body{1}.prescribed{1}.ATTR.bc='x';
febio_spec.Boundary.rigid_body{1}.prescribed{1}.ATTR.lc=1;
febio_spec.Boundary.rigid_body{1}.prescribed{1}.VAL=plateDisplacement_XYZ(1);

febio_spec.Boundary.rigid_body{1}.prescribed{2}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{1}.prescribed{2}.ATTR.lc=1;
febio_spec.Boundary.rigid_body{1}.prescribed{2}.VAL=plateDisplacement_XYZ(2);

febio_spec.Boundary.rigid_body{1}.prescribed{3}.ATTR.bc='z';
febio_spec.Boundary.rigid_body{1}.prescribed{3}.ATTR.lc=1;
febio_spec.Boundary.rigid_body{1}.prescribed{3}.VAL=plateDisplacement_XYZ(3);

febio_spec.Boundary.rigid_body{2}.ATTR.mat=4;
febio_spec.Boundary.rigid_body{2}.fixed{1}.ATTR.bc='Rx';
febio_spec.Boundary.rigid_body{2}.fixed{2}.ATTR.bc='Ry';
febio_spec.Boundary.rigid_body{2}.fixed{3}.ATTR.bc='Rz';

febio_spec.Boundary.rigid_body{2}.prescribed{1}.ATTR.bc='x';
febio_spec.Boundary.rigid_body{2}.prescribed{1}.ATTR.lc=1;
febio_spec.Boundary.rigid_body{2}.prescribed{1}.VAL=-plateDisplacement_XYZ(1);

febio_spec.Boundary.rigid_body{2}.prescribed{2}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{2}.prescribed{2}.ATTR.lc=1;
febio_spec.Boundary.rigid_body{2}.prescribed{2}.VAL=-plateDisplacement_XYZ(2);

febio_spec.Boundary.rigid_body{2}.prescribed{3}.ATTR.bc='z';
febio_spec.Boundary.rigid_body{2}.prescribed{3}.ATTR.lc=1;
febio_spec.Boundary.rigid_body{2}.prescribed{3}.VAL=-plateDisplacement_XYZ(3);

%Contact section
for q=1:1:2
    switch contactType
        case 'sticky'
            febio_spec.Contact.contact{q}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{q}.ATTR.name;
            febio_spec.Contact.contact{q}.ATTR.type='sticky';
            febio_spec.Contact.contact{q}.penalty=contactPenalty;
            febio_spec.Contact.contact{q}.laugon=0;
            febio_spec.Contact.contact{q}.tolerance=0.1;
            febio_spec.Contact.contact{q}.minaug=0;
            febio_spec.Contact.contact{q}.maxaug=10;
            febio_spec.Contact.contact{q}.snap_tol=0;
            febio_spec.Contact.contact{q}.max_traction=0;
            febio_spec.Contact.contact{q}.search_tolerance=0.1;
        case 'facet-to-facet sliding'
            febio_spec.Contact.contact{q}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{q}.ATTR.name;
            febio_spec.Contact.contact{q}.ATTR.type='facet-to-facet sliding';
            febio_spec.Contact.contact{q}.penalty=contactPenalty;
            febio_spec.Contact.contact{q}.auto_penalty=1;
            febio_spec.Contact.contact{q}.two_pass=0;
            febio_spec.Contact.contact{q}.laugon=0;
            febio_spec.Contact.contact{q}.tolerance=0.1;
            febio_spec.Contact.contact{q}.gaptol=0;
            febio_spec.Contact.contact{q}.minaug=0;
            febio_spec.Contact.contact{q}.maxaug=10;
            febio_spec.Contact.contact{q}.search_tol=0.01;
            febio_spec.Contact.contact{q}.search_radius=mean(pointSpacing)/2;
        case 'sliding_with_gaps'
            febio_spec.Contact.contact{q}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{q}.ATTR.name;
            febio_spec.Contact.contact{q}.ATTR.type='sliding_with_gaps';
            febio_spec.Contact.contact{q}.penalty=contactPenalty;
            febio_spec.Contact.contact{q}.auto_penalty=1;
            febio_spec.Contact.contact{q}.two_pass=0;
            febio_spec.Contact.contact{q}.laugon=0;
            febio_spec.Contact.contact{q}.tolerance=0.1;
            febio_spec.Contact.contact{q}.gaptol=0;
            febio_spec.Contact.contact{q}.minaug=0;
            febio_spec.Contact.contact{q}.maxaug=10;
            febio_spec.Contact.contact{q}.fric_coeff=fric_coeff;
            febio_spec.Contact.contact{q}.fric_penalty=0;
            febio_spec.Contact.contact{q}.ktmult=1;
            febio_spec.Contact.contact{q}.seg_up=0;
            febio_spec.Contact.contact{q}.search_tol=0.01;
        case 'sliding2'
            febio_spec.Contact.contact{q}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{q}.ATTR.name;
            febio_spec.Contact.contact{q}.ATTR.type='sliding2';
            febio_spec.Contact.contact{q}.penalty=contactPenalty;
            febio_spec.Contact.contact{q}.auto_penalty=1;
            febio_spec.Contact.contact{q}.two_pass=0;
            febio_spec.Contact.contact{q}.laugon=0;
            febio_spec.Contact.contact{q}.tolerance=0.1;
            febio_spec.Contact.contact{q}.gaptol=0;
            febio_spec.Contact.contact{q}.symmetric_stiffness=0;
            febio_spec.Contact.contact{q}.search_tol=0.01;
            febio_spec.Contact.contact{q}.search_radius=mean(pointSpacing)/2;
        case 'sliding-elastic'
            febio_spec.Contact.contact{q}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{q}.ATTR.name;
            febio_spec.Contact.contact{q}.ATTR.type='sliding-elastic';
            febio_spec.Contact.contact{q}.two_pass=1;
            febio_spec.Contact.contact{q}.laugon=laugon;
            febio_spec.Contact.contact{q}.tolerance=0.2;
            febio_spec.Contact.contact{q}.gaptol=0;
            febio_spec.Contact.contact{q}.minaug=minaug;
            febio_spec.Contact.contact{q}.maxaug=maxaug;
            febio_spec.Contact.contact{q}.search_tol=0.01;
            febio_spec.Contact.contact{q}.search_radius=0.1;
            febio_spec.Contact.contact{q}.symmetric_stiffness=0;
            febio_spec.Contact.contact{q}.auto_penalty=1;
            febio_spec.Contact.contact{q}.penalty=contactPenalty;
            febio_spec.Contact.contact{q}.fric_coeff=fric_coeff;
    end
end

%Output section
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{2}.VAL=1:size(V,1);

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
febioAnalysis.runMode='external';%'internal';
febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results

if runFlag==1 %i.e. a succesful run
    
    % Importing nodal displacements from a log file
    [time_mat, N_disp_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp)); %Nodal displacements
    time_mat=[0; time_mat(:)]; %Time
    
    N_disp_mat=N_disp_mat(:,2:end,:);
    sizImport=size(N_disp_mat);
    sizImport(3)=sizImport(3)+1;
    N_disp_mat_n=zeros(sizImport);
    N_disp_mat_n(:,:,2:end)=N_disp_mat;
    N_disp_mat=N_disp_mat_n;
    DN=N_disp_mat(:,:,end);
    DN_magnitude=sqrt(sum(DN(:,3).^2,2));
    V_def=V+DN;
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
    X_DEF=V_DEF(:,1,:);
    Y_DEF=V_DEF(:,2,:);
    Z_DEF=V_DEF(:,3,:);
        
    %%
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1=gpatch(Fb,V_def,DN_magnitude,'k',1); %Add graphics object to animate
    hp1.FaceColor='Interp';
    hp2=gpatch(Fp1,V_def,'kw','none',0.5); %Add graphics object to animate
    hp3=gpatch(Fp2,V_def,'kw','none',0.5); %Add graphics object to animate
    %     gpatch(Fb,V,0.5*ones(1,3),'none',0.25); %A static graphics object
    
    axisGeom(gca,fontSize);
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]);
    axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);
    camlight headlight;
    
    % Set up animation features
    animStruct.Time=time_mat; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments
        DN=N_disp_mat(:,:,qt); %Current displacement
        DN_magnitude=sqrt(sum(DN.^2,2)); %Current displacement magnitude
        V_def=V+DN; %Current nodal coordinates
                
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1 hp1 hp2 hp3]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,DN_magnitude,V_def,V_def}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct); %Initiate animation feature
    gdrawnow;
    
    
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
