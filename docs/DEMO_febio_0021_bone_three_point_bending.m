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

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting force
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

%
nSub=0;
nBlur=1;
numMaterials=25;
volumeFactor=2;

runMode='external';

%Material parameter set
c1=1e-3; %Shear-modulus-like parameter
m1=8; %Material parameter setting degree of non-linearity
k_factor=1e2; %Bulk modulus factor
k=c1*k_factor; %Bulk modulus

% FEA control settings
numTimeSteps=25; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=10; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

%Contact parameters
contactInitialOffset=0.1;
contactAlg=2;
contactPenalty=10;
switch contactAlg
    case 1
        contactType='sticky';
    case 2
        contactType='facet-to-facet sliding';
    case 3
        contactType='sliding_with_gaps';
    case 4
        contactType='sliding2';
end

zDisp=-3;

%% Prepare bone geometry

[Fs,Vs]=graphicsModels(5); %Get femur model
Vs=Vs*1000; %Scale to mm
[Fs,Vs]=triSurfRemoveThreeConnect(Fs,Vs); %remove 3 edge connected nodes and triangles

%Reorient
V_mean=mean(Vs,1);
Vs=Vs-V_mean(ones(size(Vs,1),1),:); %Center around origin
[R]=pointSetPrincipalDir(Vs); %Get rotation matrix
Vs=Vs*R; %Rotate

%Refine
[Fs,Vs]=subtri(Fs,Vs,nSub);

%Smooth
cPar.Method='HC';
cPar.n=25;
V_ini=Vs;
F_ini=Fs;
[Vs]=patchSmooth(Fs,Vs,[],cPar);

[D]=patchEdgeLengths(Fs,Vs);
voxelSize=mean(D);

%%
% Plotting surface geometry

cFigure;
subplot(1,2,1);
hold on;
gpatch(Fs,Vs,'gw','k',1);
axisGeom;
camlight headlight;

subplot(1,2,2);
hold on;
hl(1)=gpatch(F_ini,V_ini,'g','none',0.5);
hl(2)=gpatch(Fs,Vs,'r','none',0.5);
axisGeom;
camlight headlight;
legend(hl,{'Original','Smoothed'});

drawnow;

%% Create and position cylinder geometry

pointSpacing=mean(D)/2;

inputStruct.cylRadius=10;
inputStruct.numRadial=round((2*pi*inputStruct.cylRadius)./pointSpacing);
inputStruct.cylHeight=max(Vs(:,2))-min(Vs(:,2));
nh=round(inputStruct.cylHeight./pointSpacing);
nh=nh+double(iseven(nh));
inputStruct.numHeight=nh;
inputStruct.meshType='tri';
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

logicSelect=min(Vc1(:,1))<Vs(:,1) & max(Vc1(:,1))>Vs(:,1);
zOffset=max(Vs(logicSelect,3));
Vc1(:,3)=Vc1(:,3)-min(Vc1(:,3))+zOffset;

logicSelect=min(Vc2(:,1))<Vs(:,1) & max(Vc2(:,1))>Vs(:,1);
zOffset=min(Vs(logicSelect,3));
Vc2(:,3)=Vc2(:,3)-max(Vc2(:,3))+zOffset;

logicSelect=min(Vc3(:,1))<Vs(:,1) & max(Vc3(:,1))>Vs(:,1);
zOffset=min(Vs(logicSelect,3));
Vc3(:,3)=Vc3(:,3)-max(Vc3(:,3))+zOffset;

%%
% Plotting surface geometry

cFigure;
hold on;
gpatch(Fs,Vs,'kw','k',1);
gpatch(Fc,Vc1,'gw','g',1);
gpatch(Fc,Vc2,'rw','r',1);
gpatch(Fc,Vc3,'bw','b',1);

axisGeom;
camlight headlight;
drawnow

%%

% Defining the full set of possible control parameters
imOrigin=min(Vs,[],1)-voxelSize;
imMax=max(Vs,[],1)+voxelSize;
imSiz=round((imMax-imOrigin)/voxelSize);
imSiz=imSiz([2 1 3]); %Image size (x, y corresponds to j,i in image coordinates, hence the permutation)

% Using |triSurf2Im| function to convert patch data to image data
[M,~]=triSurf2Im(Fs,Vs,voxelSize,imOrigin,imSiz);

%%
% Create inner surface

L_model=(M==2); %Interior&Boundary choosen here

%Defining erosion/dilation kernel
kk=3;
p=kk-round(kk./2);
hb=zeros(3,3);
hb(2,2,2)=1;
hb(2,2,1)=1;
hb(2,2,3)=1;
hb(1,2,2)=1;
hb(3,2,2)=1;
hb(2,3,2)=1;
hb(2,1,2)=1;

for q=1:1:nBlur
    L_model_rep=zeros(size(L_model)+(2.*p));
    L_model_rep(p+(1:size(L_model,1)),p+(1:size(L_model,2)),p+(1:size(L_model,3)))=L_model;
    L_model_blur = convn(double(L_model_rep),hb,'valid');
    L_model=L_model_blur>=(sum(hb(:)));
end

[Fs2,Vs2,~]=im2patch(L_model,L_model,'vb',voxelSize*ones(1,3));
Fs2=[Fs2(:,[1 2 3]); Fs2(:,[3 4 1])];
Vs2=Vs2+imOrigin(ones(size(Vs2,1),1),:);

cPar.Method='HC';
cPar.n=25;
[Vs2]=patchSmooth(Fs2,Vs2,[],cPar);

%%
% Visualize
cFigure;
hold on;
gpatch(Fs,Vs,'gw','none',0.5);
gpatch(Fs2,Vs2,'bw','k',1);
axisGeom;
camlight headlight;
drawnow;

%% Get interior point

[M2,~]=triSurf2Im(Fs2,Vs2,voxelSize,imOrigin,imSiz);
L=M>0 & M2==0;

[indInternal]=getInnerVoxel(L,6,0);

[ii,jj,kk]=ind2sub(size(L),indInternal);
V_in1=nan(1,3);
[V_in1(:,1),V_in1(:,2),V_in1(:,3)]=im2cart(ii,jj,kk,voxelSize);
V_in1=V_in1+imOrigin;

[V_in2]=getInnerPoint(Fs2,Vs2,[],[],1);

%%
% Visualize
cFigure;
hold on;
gpatch(Fs,Vs,'gw','none',0.5);
gpatch(Fs2,Vs2,'bw','none',0.5);
plotV(V_in1,'k.','MarkerSize',25);
axisGeom;
camlight headlight;
drawnow;

%%

[F,V,C]=joinElementSets({Fs,Fs2},{Vs,Vs2});

%%
% Visualize
cFigure;
hold on;
gpatch(F,V,C,'k',0.5);
colormap(gjet(2)); icolorbar;
axisGeom;
camlight headlight;
drawnow;

%%
% CREATING THE INPUT STRUCTURE

[regionA]=tetVolMeanEst(F,V); %Volume estimate for regular tets

stringOpt='-pq1.2AaY';

inputStruct.stringOpt=stringOpt;
inputStruct.Faces=F;
inputStruct.Nodes=V;
inputStruct.holePoints=V_in2;
inputStruct.faceBoundaryMarker=C; %Face boundary markers
inputStruct.regionPoints=V_in1; %region points
inputStruct.regionA=regionA*volumeFactor;
inputStruct.minRegionMarker=2; %Minimum region marker

% Mesh model using tetrahedral elements using tetGen
[meshOutput]=runTetGen(inputStruct); %Run tetGen

% Access model element and patch data
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
V=meshOutput.nodes;
CE=meshOutput.elementMaterialID;
E=meshOutput.elements;

%% Visualizing mesh using |meshView|, see also |anim8|

meshView(meshOutput);

%%

% D=minDist(V,Vs);
% D=vertexToFaceMeasure(E,D);
% D=D-min(D(:));
% D=D./max(D(:));
% D=abs(D-1);
% D=round(1+(D*(numMaterials-1)));
%
% meshOutput2=meshOutput;
% meshOutput2.elementMaterialID=D;
%
% %%
% meshView(meshOutput2);

%% Joining node sets

[E_F,V,C]=joinElementSets({E,Fc,Fc,Fc},{V,Vc1,Vc2,Vc3});
% E=E_F{1};
% Fb
Fc1=E_F{2};
Fc2=E_F{3};
Fc3=E_F{4};

%%
% Visualize

Fb_slave=Fb(Cb==1,:);

cFigure;
hold on;

gpatch(Fb_slave,V,'kw','k',0.5);
gpatch(Fc1,V,'rw','r',1);
gpatch(Fc2,V,'gw','g',1);
gpatch(Fc3,V,'bw','b',1);

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
febio_spec.Control.title='Cube analysis';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax;
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;
febio_spec.Control.max_refs=max_refs;
febio_spec.Control.max_ups=max_ups;

%Material section
febio_spec.Material.material{1}.ATTR.type='Ogden';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c1=c1;
febio_spec.Material.material{1}.m1=m1;
febio_spec.Material.material{1}.c2=c1;
febio_spec.Material.material{1}.m2=-m1;
febio_spec.Material.material{1}.k=k;

febio_spec.Material.material{2}.ATTR.type='rigid body';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.density=1;
febio_spec.Material.material{2}.center_of_mass=mean(Vc1,1);

febio_spec.Material.material{3}.ATTR.type='rigid body';
febio_spec.Material.material{3}.ATTR.id=3;
febio_spec.Material.material{3}.density=1;
febio_spec.Material.material{3}.center_of_mass=mean(Vc2,1);

febio_spec.Material.material{4}.ATTR.type='rigid body';
febio_spec.Material.material{4}.ATTR.id=4;
febio_spec.Material.material{4}.density=1;
febio_spec.Material.material{4}.center_of_mass=mean(Vc3,1);

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='tet4'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set
febio_spec.Geometry.Elements{1}.ATTR.name='bone'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=E;

febio_spec.Geometry.Elements{2}.ATTR.type='tri3'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat=2; %material index for this set
febio_spec.Geometry.Elements{2}.ATTR.name='bar1'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id=size(E,1)+(1:1:size(Fc1,1))'; %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL=Fc1;

febio_spec.Geometry.Elements{3}.ATTR.type='tri3'; %Element type of this set
febio_spec.Geometry.Elements{3}.ATTR.mat=3; %material index for this set
febio_spec.Geometry.Elements{3}.ATTR.name='bar2'; %Name of the element set
febio_spec.Geometry.Elements{3}.elem.ATTR.id=size(E,1)+(1:1:size(Fc2,1))'; %Element id's
febio_spec.Geometry.Elements{3}.elem.VAL=Fc2;

febio_spec.Geometry.Elements{4}.ATTR.type='tri3'; %Element type of this set
febio_spec.Geometry.Elements{4}.ATTR.mat=4; %material index for this set
febio_spec.Geometry.Elements{4}.ATTR.name='bar3'; %Name of the element set
febio_spec.Geometry.Elements{4}.elem.ATTR.id=size(E,1)+(1:1:size(Fc3,1))'; %Element id's
febio_spec.Geometry.Elements{4}.elem.VAL=Fc3;

% -> NodeSets
% febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList';
% febio_spec.Geometry.NodeSet{1}.VAL=bcSupportList(:);

% -> Surfaces
febio_spec.Geometry.Surface{1}.ATTR.name='contact_master1';
febio_spec.Geometry.Surface{1}.tri3.ATTR.lid=(1:1:size(Fc1,1))';
febio_spec.Geometry.Surface{1}.tri3.VAL=Fc1;

febio_spec.Geometry.Surface{2}.ATTR.name='contact_master2';
febio_spec.Geometry.Surface{2}.tri3.ATTR.lid=(1:1:size(Fc2,1))';
febio_spec.Geometry.Surface{2}.tri3.VAL=Fc2;

febio_spec.Geometry.Surface{3}.ATTR.name='contact_master3';
febio_spec.Geometry.Surface{3}.tri3.ATTR.lid=(1:1:size(Fc3,1))';
febio_spec.Geometry.Surface{3}.tri3.VAL=Fc3;

febio_spec.Geometry.Surface{4}.ATTR.name='contact_slave';
febio_spec.Geometry.Surface{4}.tri3.ATTR.lid=(1:1:size(Fb_slave,1))';
febio_spec.Geometry.Surface{4}.tri3.VAL=Fb_slave;

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
% -> Fix boundary conditions
% febio_spec.Boundary.fix{1}.ATTR.bc='x';
% febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
% febio_spec.Boundary.fix{2}.ATTR.bc='y';
% febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
% febio_spec.Boundary.fix{3}.ATTR.bc='z';
% febio_spec.Boundary.fix{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;

% -> Prescribed boundary conditions on the rigid body
febio_spec.Boundary.rigid_body{1}.ATTR.mat=2;
febio_spec.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='x';
febio_spec.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='Rx';
febio_spec.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Ry';
febio_spec.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Rz';
febio_spec.Boundary.rigid_body{1}.prescribed.ATTR.bc='z';
febio_spec.Boundary.rigid_body{1}.prescribed.ATTR.lc=1;
febio_spec.Boundary.rigid_body{1}.prescribed.VAL=zDisp;

febio_spec.Boundary.rigid_body{2}.ATTR.mat=3;
febio_spec.Boundary.rigid_body{2}.fixed{1}.ATTR.bc='x';
febio_spec.Boundary.rigid_body{2}.fixed{2}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{2}.fixed{3}.ATTR.bc='z';
febio_spec.Boundary.rigid_body{2}.fixed{4}.ATTR.bc='Rx';
febio_spec.Boundary.rigid_body{2}.fixed{5}.ATTR.bc='Ry';
febio_spec.Boundary.rigid_body{2}.fixed{6}.ATTR.bc='Rz';

febio_spec.Boundary.rigid_body{3}.ATTR.mat=4;
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
            febio_spec.Contact.contact{qc}.tolerance=0.1;
            febio_spec.Contact.contact{qc}.minaug=0;
            febio_spec.Contact.contact{qc}.maxaug=10;
            febio_spec.Contact.contact{qc}.snap_tol=0;
            febio_spec.Contact.contact{qc}.max_traction=0;
            febio_spec.Contact.contact{qc}.search_tolerance=0.1;
        case 'facet-to-facet sliding'
            febio_spec.Contact.contact{qc}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{qc}.ATTR.name;
            febio_spec.Contact.contact{qc}.ATTR.type='facet-to-facet sliding';
            febio_spec.Contact.contact{qc}.penalty=contactPenalty;
            febio_spec.Contact.contact{qc}.auto_penalty=1;
            febio_spec.Contact.contact{qc}.two_pass=0;
            febio_spec.Contact.contact{qc}.laugon=0;
            febio_spec.Contact.contact{qc}.tolerance=0.1;
            febio_spec.Contact.contact{qc}.gaptol=0;
            febio_spec.Contact.contact{qc}.minaug=0;
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
febioAnalysis.runMode=runMode;
febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=3; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results

if runFlag==1 || runFlag==0 %i.e. a succesful run
    
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
    [CF]=vertexToFaceMeasure(Fb,DN_magnitude);
    
    %
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1=gpatch(Fb,V_def,CF,'k',1); %Add graphics object to animate
    hp2=gpatch([Fc1;Fc2;Fc3],V_def,'kw','none',0.5);
    
    gpatch(Fb,V,0.5*ones(1,3),'none',0.25); %A static graphics object
    
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
        [CF]=vertexToFaceMeasure(Fb,DN_magnitude); %Current color data to use
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1 hp1 hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,CF,V_def}; %Property values for to set in order to animate
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
% Copyright (C) 2018  Kevin Mattheus Moerman
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
