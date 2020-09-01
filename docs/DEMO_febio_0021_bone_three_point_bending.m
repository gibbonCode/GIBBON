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
% * febio_spec version 2.5
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
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress

%Geometric parameters
corticalThickness=3; %Thickness used for cortical material definition
volumeFactor=10; %Factor to scale desired volume for interior elements w.r.t. boundary elements

%Material parameter set
E_youngs1=17000; %Youngs modulus
nu1=0.25; %Poissons ratio

% Cancellous bone
E_youngs2=1500; %Youngs modulus
nu2=0.25; %Poissons ratio

% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=15; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
min_residual=1e-20;
symmetric_stiffness=1;
runMode='internal';

%Contact parameters
contactInitialOffset=0.01;
contactAlg=1;
contactPenalty=1000;
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

zDisp=-6-2*contactInitialOffset;

%% Prepare bone geometry

% Import bone surface model
[stlStruct] = import_STL(fullfile(pathNameSTL,'femur_iso.stl'));
F_bone=stlStruct.solidFaces{1}; %Faces
V_bone=stlStruct.solidVertices{1}; %Vertices
V_bone=V_bone*1000; %Scale to mm
[F_bone,V_bone]=mergeVertices(F_bone,V_bone); % Merging nodes
F_bone=fliplr(F_bone);

%Reorient
V_mean=mean(V_bone,1);
V_bone=V_bone-V_mean(ones(size(V_bone,1),1),:); %Center around origin
[R]=pointSetPrincipalDir(V_bone); %Get rotation matrix
V_bone=V_bone*R; %Rotate

[D]=patchEdgeLengths(F_bone,V_bone);
voxelSize=mean(D);

%% Create and position cylinder geometry

pointSpacing=mean(D)/2;

inputStruct.cylRadius=10;
inputStruct.numRadial=round((2*pi*inputStruct.cylRadius)./pointSpacing);
inputStruct.cylHeight=max(V_bone(:,2))-min(V_bone(:,2));
nh=round(inputStruct.cylHeight./pointSpacing);
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

%% Visualize bone surface

cFigure; hold on;
gpatch(F_bone,V_bone,'w','k',1);
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

%% Define bone contact surface
Fb_slave=fliplr(Fb);

%%
% Visualize

cFigure; hold on;

gpatch(Fb_slave,V,'w','k',1);
patchNormPlot(Fb_slave,V);
gpatch(Fc1,V,'rw','r',1);
patchNormPlot(Fc1,V);
gpatch(Fc2,V,'gw','g',1);
patchNormPlot(Fc2,V);
gpatch(Fc3,V,'bw','b',1);
patchNormPlot(Fc3,V);
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
febio_spec.Control.symmetric_stiffness=symmetric_stiffness; 
febio_spec.Control.min_residual=min_residual;

%Material section
febio_spec.Material.material{1}.ATTR.type='neo-Hookean';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.E=E_youngs1;
febio_spec.Material.material{1}.v=nu1;

febio_spec.Material.material{2}.ATTR.type='neo-Hookean';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.E=E_youngs2;
febio_spec.Material.material{2}.v=nu2;

% febio_spec.Material.material{1}.ATTR.type='Ogden';
% febio_spec.Material.material{1}.ATTR.id=1;
% febio_spec.Material.material{1}.c1=1e-3;
% febio_spec.Material.material{1}.m1=2;
% febio_spec.Material.material{1}.c2=1e-3;
% febio_spec.Material.material{1}.m2=2;
% febio_spec.Material.material{1}.k=0.1;
% 
% febio_spec.Material.material{2}.ATTR.type='Ogden';
% febio_spec.Material.material{2}.ATTR.id=2;
% febio_spec.Material.material{2}.c1=1e-3;
% febio_spec.Material.material{2}.m1=2;
% febio_spec.Material.material{2}.c2=1e-3;
% febio_spec.Material.material{2}.m2=2;
% febio_spec.Material.material{2}.k=0.1;

febio_spec.Material.material{3}.ATTR.type='rigid body';
febio_spec.Material.material{3}.ATTR.id=3;
febio_spec.Material.material{3}.density=1;
febio_spec.Material.material{3}.center_of_mass=mean(Vc1,1);

febio_spec.Material.material{4}.ATTR.type='rigid body';
febio_spec.Material.material{4}.ATTR.id=4;
febio_spec.Material.material{4}.density=1;
febio_spec.Material.material{4}.center_of_mass=mean(Vc2,1);

febio_spec.Material.material{5}.ATTR.type='rigid body';
febio_spec.Material.material{5}.ATTR.id=5;
febio_spec.Material.material{5}.density=1;
febio_spec.Material.material{5}.center_of_mass=mean(Vc3,1);

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

febio_spec.Geometry.Elements{3}.ATTR.type='quad4'; %Element type of this set
febio_spec.Geometry.Elements{3}.ATTR.mat=3; %material index for this set
febio_spec.Geometry.Elements{3}.ATTR.name='bar1'; %Name of the element set
febio_spec.Geometry.Elements{3}.elem.ATTR.id=size(E1,1)+size(E2,1)+(1:1:size(Fc1,1))'; %Element id's
febio_spec.Geometry.Elements{3}.elem.VAL=Fc1;

febio_spec.Geometry.Elements{4}.ATTR.type='quad4'; %Element type of this set
febio_spec.Geometry.Elements{4}.ATTR.mat=4; %material index for this set
febio_spec.Geometry.Elements{4}.ATTR.name='bar2'; %Name of the element set
febio_spec.Geometry.Elements{4}.elem.ATTR.id=size(E1,1)+size(E2,1)+size(Fc1,1)+(1:1:size(Fc2,1))'; %Element id's
febio_spec.Geometry.Elements{4}.elem.VAL=Fc2;

febio_spec.Geometry.Elements{5}.ATTR.type='quad4'; %Element type of this set
febio_spec.Geometry.Elements{5}.ATTR.mat=5; %material index for this set
febio_spec.Geometry.Elements{5}.ATTR.name='bar3'; %Name of the element set
febio_spec.Geometry.Elements{5}.elem.ATTR.id=size(E1,1)+size(E2,1)+size(Fc1,1)+size(Fc2,1)+(1:1:size(Fc3,1))'; %Element id's
febio_spec.Geometry.Elements{5}.elem.VAL=Fc3;

% -> Surfaces
febio_spec.Geometry.Surface{1}.ATTR.name='contact_master1';
febio_spec.Geometry.Surface{1}.quad4.ATTR.lid=(1:1:size(Fc1,1))';
febio_spec.Geometry.Surface{1}.quad4.VAL=Fc1;

febio_spec.Geometry.Surface{2}.ATTR.name='contact_master2';
febio_spec.Geometry.Surface{2}.quad4.ATTR.lid=(1:1:size(Fc2,1))';
febio_spec.Geometry.Surface{2}.quad4.VAL=Fc2;

febio_spec.Geometry.Surface{3}.ATTR.name='contact_master3';
febio_spec.Geometry.Surface{3}.quad4.ATTR.lid=(1:1:size(Fc3,1))';
febio_spec.Geometry.Surface{3}.quad4.VAL=Fc3;

febio_spec.Geometry.Surface{4}.ATTR.name='contact_slave';
febio_spec.Geometry.Surface{4}.quad4.ATTR.lid=(1:1:size(Fb_slave,1))';
febio_spec.Geometry.Surface{4}.quad4.VAL=Fb_slave;

% -> Surface pairs
febio_spec.Geometry.SurfacePair{1}.ATTR.name='Contact1';
febio_spec.Geometry.SurfacePair{1}.master.ATTR.surface=febio_spec.Geometry.Surface{1}.ATTR.name;
febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface=febio_spec.Geometry.Surface{4}.ATTR.name;

febio_spec.Geometry.SurfacePair{2}.ATTR.name='Contact2';
febio_spec.Geometry.SurfacePair{2}.master.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;
febio_spec.Geometry.SurfacePair{2}.slave.ATTR.surface=febio_spec.Geometry.Surface{4}.ATTR.name;

febio_spec.Geometry.SurfacePair{3}.ATTR.name='Contact3';
febio_spec.Geometry.SurfacePair{3}.master.ATTR.surface=febio_spec.Geometry.Surface{3}.ATTR.name;
febio_spec.Geometry.SurfacePair{3}.slave.ATTR.surface=febio_spec.Geometry.Surface{4}.ATTR.name;

%Boundary condition section

% -> Prescribed boundary conditions on the rigid body
febio_spec.Boundary.rigid_body{1}.ATTR.mat=3;
febio_spec.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='x';
febio_spec.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='Rx';
febio_spec.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Ry';
febio_spec.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Rz';
febio_spec.Boundary.rigid_body{1}.prescribed.ATTR.bc='z';
febio_spec.Boundary.rigid_body{1}.prescribed.ATTR.lc=1;
febio_spec.Boundary.rigid_body{1}.prescribed.VAL=zDisp;

febio_spec.Boundary.rigid_body{2}.ATTR.mat=4;
febio_spec.Boundary.rigid_body{2}.fixed{1}.ATTR.bc='x';
febio_spec.Boundary.rigid_body{2}.fixed{2}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{2}.fixed{3}.ATTR.bc='z';
febio_spec.Boundary.rigid_body{2}.fixed{4}.ATTR.bc='Rx';
febio_spec.Boundary.rigid_body{2}.fixed{5}.ATTR.bc='Ry';
febio_spec.Boundary.rigid_body{2}.fixed{6}.ATTR.bc='Rz';

febio_spec.Boundary.rigid_body{3}.ATTR.mat=5;
febio_spec.Boundary.rigid_body{3}.fixed{1}.ATTR.bc='x';
febio_spec.Boundary.rigid_body{3}.fixed{2}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{3}.fixed{3}.ATTR.bc='z';
febio_spec.Boundary.rigid_body{3}.fixed{4}.ATTR.bc='Rx';
febio_spec.Boundary.rigid_body{3}.fixed{5}.ATTR.bc='Ry';
febio_spec.Boundary.rigid_body{3}.fixed{6}.ATTR.bc='Rz';

for qc=1:1:3
    
    %Contact section
    switch contactType
        case 'sticky'
            febio_spec.Contact.contact{qc}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{qc}.ATTR.name;
            febio_spec.Contact.contact{qc}.ATTR.type='sticky';
            febio_spec.Contact.contact{qc}.penalty=contactPenalty;
            febio_spec.Contact.contact{qc}.laugon=0;
            febio_spec.Contact.contact{qc}.tolerance=0.2;
            febio_spec.Contact.contact{qc}.minaug=1;
            febio_spec.Contact.contact{qc}.maxaug=10;
            febio_spec.Contact.contact{qc}.snap_tol=0;
            febio_spec.Contact.contact{qc}.max_traction=0;
            febio_spec.Contact.contact{qc}.search_tolerance=0.01;
        case 'facet-to-facet sliding'
            febio_spec.Contact.contact{qc}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{qc}.ATTR.name;
            febio_spec.Contact.contact{qc}.ATTR.type='facet-to-facet sliding';
            febio_spec.Contact.contact{qc}.penalty=contactPenalty;
            febio_spec.Contact.contact{qc}.auto_penalty=1;
            febio_spec.Contact.contact{qc}.two_pass=1;
            febio_spec.Contact.contact{qc}.laugon=1;
            febio_spec.Contact.contact{qc}.tolerance=0.1;
            febio_spec.Contact.contact{qc}.gaptol=0;
            febio_spec.Contact.contact{qc}.minaug=1;
            febio_spec.Contact.contact{qc}.maxaug=10;
            febio_spec.Contact.contact{qc}.search_tol=0.01;
            febio_spec.Contact.contact{qc}.search_radius=mean(pointSpacing)/2;
        case 'sliding_with_gaps'
            febio_spec.Contact.contact{qc}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{qc}.ATTR.name;
            febio_spec.Contact.contact{qc}.ATTR.type='sliding_with_gaps';
            febio_spec.Contact.contact{qc}.penalty=contactPenalty;
            febio_spec.Contact.contact{qc}.auto_penalty=1;
            febio_spec.Contact.contact{qc}.two_pass=0;
            febio_spec.Contact.contact{qc}.laugon=0;
            febio_spec.Contact.contact{qc}.tolerance=0.1;
            febio_spec.Contact.contact{qc}.gaptol=0;
            febio_spec.Contact.contact{qc}.minaug=0;
            febio_spec.Contact.contact{qc}.maxaug=10;
            febio_spec.Contact.contact{qc}.fric_coeff=0;
            febio_spec.Contact.contact{qc}.fric_penalty=0;
            febio_spec.Contact.contact{qc}.ktmult=1;
            febio_spec.Contact.contact{qc}.seg_up=0;
            febio_spec.Contact.contact{qc}.search_tol=0.01;
        case 'sliding2'
            febio_spec.Contact.contact{qc}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{qc}.ATTR.name;
            febio_spec.Contact.contact{qc}.ATTR.type='sliding2';
            febio_spec.Contact.contact{qc}.penalty=contactPenalty;
            febio_spec.Contact.contact{qc}.auto_penalty=1;
            febio_spec.Contact.contact{qc}.two_pass=0;
            febio_spec.Contact.contact{qc}.laugon=0;
            febio_spec.Contact.contact{qc}.tolerance=0.1;
            febio_spec.Contact.contact{qc}.gaptol=0;
            febio_spec.Contact.contact{qc}.symmetric_stiffness=0;
            febio_spec.Contact.contact{qc}.search_tol=0.01;
            febio_spec.Contact.contact{qc}.search_radius=mean(pointSpacings)/2;
        case 'sliding-elastic'
            febio_spec.Contact.contact{qc}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{qc}.ATTR.name;
            febio_spec.Contact.contact{qc}.ATTR.type='sliding-elastic';
            febio_spec.Contact.contact{qc}.two_pass=1;
            febio_spec.Contact.contact{qc}.laugon=1;
            febio_spec.Contact.contact{qc}.tolerance=0.2;
            febio_spec.Contact.contact{qc}.gaptol=0;
            febio_spec.Contact.contact{qc}.minaug=1;
            febio_spec.Contact.contact{qc}.maxaug=10;
            febio_spec.Contact.contact{qc}.search_tol=0.01;
            febio_spec.Contact.contact{qc}.search_radius=1;
            febio_spec.Contact.contact{qc}.symmetric_stiffness=0;
            febio_spec.Contact.contact{qc}.auto_penalty=1;
            febio_spec.Contact.contact{qc}.penalty=contactPenalty;
            febio_spec.Contact.contact{qc}.fric_coeff=0.5;
    end
end

%Output section
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

febio_spec.Output.logfile.rigid_body_data{1}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.rigid_body_data{1}.ATTR.data='Fx;Fy;Fz';
febio_spec.Output.logfile.rigid_body_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.rigid_body_data{1}.VAL=3;

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='s1';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.element_data{1}.VAL=1:size(E,1);

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
febioAnalysis.runMode=runMode;
febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results 

if runFlag==1 %i.e. a succesful run
    
    %% 
    % Importing nodal displacements from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
    
    %Access data
    N_disp_mat=dataStruct.data; %Displacement
    timeVec=dataStruct.time; %Time
    
    %Create deformed coordinate set
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
               
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    DN_magnitude=sqrt(sum(N_disp_mat(:,:,end).^2,2)); %Current displacement magnitude
        
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('Displacement magnitude [mm]','Interpreter','Latex')
    hp=gpatch(Fb,V_DEF(:,:,end),DN_magnitude,'k',1); %Add graphics object to animate
    hp.FaceColor='interp';
    
    hp2=gpatch([Fc1;Fc2;Fc3],V,'w','none',1); 
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]);    
    axis(axisLim(V_DEF)); %Set axis limits statically    
    camlight headlight;        
        
    % Set up animation features
    animStruct.Time=timeVec; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        DN_magnitude=sqrt(sum(N_disp_mat(:,:,qt).^2,2)); %Current displacement magnitude
                
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),DN_magnitude,V_DEF(:,:,qt)}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;
            
    %%
    % Importing element stress from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),1,1);
    
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
    hp2=gpatch([Fc1;Fc2;Fc3],V,'w','none',0.5); 
    
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
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_force),1,1);
    
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
