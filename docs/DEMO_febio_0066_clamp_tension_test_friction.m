%% DEMO_febio_0066_clamp_tension_test_friction
% Below is a demonstration for: 
% 1) The creation of an FEBio model for clamped tensile testing
% 2) The use of multiple steps
% 4) Running an FEBio job with MATLAB
% 5) Importing FEBio results into MATLAB

%% Keywords
%
% * febio_spec version 2.5
% * febio, FEBio
% * tension, tensile
% * displacement control, displacement boundary condition
% * Hexahedral hex8
% * static, solid
% * hyperelastic, Ogden
% * displacement logfile
% * stress logfile

%%

clear; close all; clc;

%% Plot settings
fontSize=25;

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
pointSpacingStrip=1; 
sampleWidth=10;
sampleThickness=2; 
sampleGripGripHeight=sampleWidth.*2;
sampleClampedHeight=10;

pointSpacingClamp=pointSpacingStrip/2; 
clampExtendAmount=sampleWidth/2;

tensileStrain=0.3;
clampCompressiveStrain=0.3;
initialSpacing=pointSpacingStrip/10;

%Initial material parameter set
c1=1e-3;
m1=2;
k_factor=100;
k=c1*k_factor;

% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=12; %Optimum number of iterations
max_retries=25; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

min_residual=1e-20;
symmetric_stiffness=0;
runMode='external'; %'internal' or 'external'

%Contact parameters
contactPenalty=10;
laugon=0;
minaug=1;
maxaug=10;
fric_coeff=0.5;

%% Computing derived parameters
sampleHeight=sampleGripGripHeight+2*sampleClampedHeight;
numElementsWidthStrip=ceil(sampleWidth/pointSpacingStrip);
numElementsWidthStrip=numElementsWidthStrip+iseven(numElementsWidthStrip); %Force uneven so there is a middle element
numElementsThicknessStrip=ceil(round(sampleThickness/pointSpacingStrip));
numElementsHeightStrip=ceil(round(sampleHeight/pointSpacingStrip));

clampWidth=sampleWidth+(2*clampExtendAmount);
clampHeight=sampleClampedHeight+(clampExtendAmount);
numElementsWidthClamp=ceil(clampWidth/pointSpacingClamp);
numElementsHeightClamp=ceil(clampHeight/pointSpacingClamp);

clampCompressiveDisplacement=(sampleThickness.*clampCompressiveStrain)/2+initialSpacing;
clampTensionDisplacement=(sampleGripGripHeight.*tensileStrain);

%% Creating strip region
% The region consists of three "boxes" which define the upper and lower
% clamped regions as well as the central region. 

%Create box 1
boxDim=[sampleWidth sampleThickness sampleHeight]; %Dimensions
boxElem=[numElementsWidthStrip numElementsThicknessStrip numElementsHeightStrip]; %Number of elements
[boxMesh]=hexMeshBox(boxDim,boxElem);
E=boxMesh.E;
V=boxMesh.V;
F=boxMesh.F;
Fb=boxMesh.Fb;
faceBoundaryMarker=boxMesh.faceBoundaryMarker;

% R1=euler2DCM([0 (8/180)*pi 0]);
% V=V*R1;

%%
% Plotting surface models
cFigure; hold on;
title('The strip mesh','FontSize',fontSize);
gpatch(Fb,V,faceBoundaryMarker);
axisGeom(gca,fontSize);
colormap(gca,gjet(250)); icolorbar; 
drawnow; 

%% Define rigid clamping surfaces

[Fc,Vc]=quadPlate([clampWidth clampHeight],[numElementsWidthClamp numElementsHeightClamp]);
R=euler2DCM([-pi/2 0 0]);
Vc=Vc*R;


Fc1=Fc;
Vc1=Vc;
Vc1(:,2)=Vc1(:,2)+sampleThickness/2+initialSpacing;
Vc1(:,3)=Vc1(:,3)-min(Vc1(:,3))+sampleGripGripHeight/2;

Vc2=Vc1;
Vc2(:,2)=-Vc2(:,2);
Fc2=fliplr(Fc);

Fc3=fliplr(Fc2);
Vc3=-Vc2;

Fc4=fliplr(Fc1);
Vc4=-Vc1;

%% 
% Visualize clamping surfaces 

cFigure; hold on;
title('Clamping surfaces','FontSize',fontSize);
gpatch(Fb,V,'kw','none',0.25);
gpatch(Fc1,Vc1,'rw','k',1);
patchNormPlot(Fc1,Vc1);
gpatch(Fc2,Vc2,'gw','k',1);
patchNormPlot(Fc2,Vc2);
gpatch(Fc3,Vc3,'bw','k',1);
patchNormPlot(Fc3,Vc3);
gpatch(Fc4,Vc4,'yw','k',1);
patchNormPlot(Fc4,Vc4);
axisGeom(gca,fontSize);
camlight headlight;
drawnow;  

%% Join node sets

Fc1=Fc1+size(V,1);
Fc2=Fc2+size(V,1)+1*size(Vc1,1);
Fc3=Fc3+size(V,1)+2*size(Vc1,1);
Fc4=Fc4+size(V,1)+3*size(Vc1,1);

V=[V;Vc1;Vc2;Vc3;Vc4];

Fc=[Fc1;Fc2;Fc3;Fc4];
Cc=repmat(1:4,size(Fc1,1),1);
Cc=Cc(:);

%%

Fb1=Fb(faceBoundaryMarker==3,:);
Fb2=Fb(faceBoundaryMarker==4,:);

%%
% Visualize boundary conditions

cFigure; hold on;
title('Complete model','FontSize',fontSize);

gpatch(Fb,V,'kw','none',0.25);

gpatch(Fb1,V,'cw','k',1);
gpatch(Fb2,V,'kw','k',1);
patchNormPlot(Fb1,V);
patchNormPlot(Fb2,V);

gpatch(Fc1,V,'rw','k',1);
gpatch(Fc2,V,'gw','k',1);
gpatch(Fc3,V,'bw','k',1);
gpatch(Fc4,V,'yw','k',1);

axisGeom(gca,fontSize);
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

%Create control structure for use by all steps
stepStruct.Control.analysis.ATTR.type='static';
stepStruct.Control.time_steps=numTimeSteps;
stepStruct.Control.step_size=1/numTimeSteps;
stepStruct.Control.time_stepper.dtmin=dtmin;
stepStruct.Control.time_stepper.dtmax=dtmax; 
stepStruct.Control.time_stepper.max_retries=max_retries;
stepStruct.Control.time_stepper.opt_iter=opt_iter;
stepStruct.Control.max_refs=max_refs;
stepStruct.Control.max_ups=max_ups;
stepStruct.Control.symmetric_stiffness=symmetric_stiffness;
stepStruct.Control.min_residual=min_residual;

%Add template based default settings to proposed control section
[stepStruct.Control]=structComplete(stepStruct.Control,febio_spec.Control,1); %Complement provided with default if missing

%Remove control field (part of template) since step specific control sections are used
febio_spec=rmfield(febio_spec,'Control'); 

febio_spec.Step{1}.ATTR.id=1;
febio_spec.Step{1}.Control=stepStruct.Control;
febio_spec.Step{2}.ATTR.id=2;
febio_spec.Step{2}.Control=stepStruct.Control;

%Material section
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.ATTR.name='Normal material';
febio_spec.Material.material{1}.ATTR.type='Ogden';
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

febio_spec.Material.material{5}.ATTR.type='rigid body';
febio_spec.Material.material{5}.ATTR.id=5;
febio_spec.Material.material{5}.density=1;
febio_spec.Material.material{5}.center_of_mass=mean(Vc4,1);

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='hex8'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name='Strip'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.VAL=E;
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:size(E,1))'; %Element material id's

febio_spec.Geometry.Elements{2}.ATTR.type='quad4'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat=2; %material index for this set 
febio_spec.Geometry.Elements{2}.ATTR.name='Grip1'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.VAL=Fc1;
febio_spec.Geometry.Elements{2}.elem.ATTR.id=size(E,1)+(1:size(Fc1,1))'; %Element material id's

febio_spec.Geometry.Elements{3}.ATTR.type='quad4'; %Element type of this set
febio_spec.Geometry.Elements{3}.ATTR.mat=3; %material index for this set 
febio_spec.Geometry.Elements{3}.ATTR.name='Grip2'; %Name of the element set
febio_spec.Geometry.Elements{3}.elem.VAL=Fc2;
febio_spec.Geometry.Elements{3}.elem.ATTR.id=size(E,1)+size(Fc1,1)+(1:size(Fc2,1))'; %Element material id's

febio_spec.Geometry.Elements{4}.ATTR.type='quad4'; %Element type of this set
febio_spec.Geometry.Elements{4}.ATTR.mat=4; %material index for this set 
febio_spec.Geometry.Elements{4}.ATTR.name='Grip3'; %Name of the element set
febio_spec.Geometry.Elements{4}.elem.VAL=Fc3;
febio_spec.Geometry.Elements{4}.elem.ATTR.id=size(E,1)+2*size(Fc1,1)+(1:size(Fc3,1))'; %Element material id's

febio_spec.Geometry.Elements{5}.ATTR.type='quad4'; %Element type of this set
febio_spec.Geometry.Elements{5}.ATTR.mat=5; %material index for this set 
febio_spec.Geometry.Elements{5}.ATTR.name='Grip4'; %Name of the element set
febio_spec.Geometry.Elements{5}.elem.VAL=Fc4;
febio_spec.Geometry.Elements{5}.elem.ATTR.id=size(E,1)+3*size(Fc1,1)+(1:size(Fc4,1))'; %Element material id's

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

febio_spec.Geometry.Surface{4}.ATTR.name='contact_master4';
febio_spec.Geometry.Surface{4}.quad4.ATTR.lid=(1:1:size(Fc4,1))';
febio_spec.Geometry.Surface{4}.quad4.VAL=Fc4;

febio_spec.Geometry.Surface{5}.ATTR.name='contact_slave1';
febio_spec.Geometry.Surface{5}.quad4.ATTR.lid=(1:1:size(Fb1,1))';
febio_spec.Geometry.Surface{5}.quad4.VAL=Fb1;

febio_spec.Geometry.Surface{6}.ATTR.name='contact_slave2';
febio_spec.Geometry.Surface{6}.quad4.ATTR.lid=(1:1:size(Fb2,1))';
febio_spec.Geometry.Surface{6}.quad4.VAL=Fb2;

% -> Surface pairs
for q=1:1:4
    if iseven(q)
        indStrip=5;
    else
        indStrip=6;
    end    
    febio_spec.Geometry.SurfacePair{q}.ATTR.name=['Contact_',num2str(q)];
    febio_spec.Geometry.SurfacePair{q}.master.ATTR.surface=febio_spec.Geometry.Surface{q}.ATTR.name;
    febio_spec.Geometry.SurfacePair{q}.slave.ATTR.surface=febio_spec.Geometry.Surface{indStrip}.ATTR.name;    
end

%Boundary condition section 
% -> Prescribe boundary conditions
%STEP 1 Clamp compression

%Set 1
% -> Prescribed boundary conditions on the rigid body
febio_spec.Step{1}.Boundary.rigid_body{1}.ATTR.mat=2;
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='x';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='z';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='Rx';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Ry';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Rz';
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.ATTR.bc='y';
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.ATTR.lc=1;
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.VAL=-clampCompressiveDisplacement;

febio_spec.Step{1}.Boundary.rigid_body{2}.ATTR.mat=3;
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{1}.ATTR.bc='x';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{2}.ATTR.bc='z';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{3}.ATTR.bc='Rx';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{4}.ATTR.bc='Ry';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{5}.ATTR.bc='Rz';
febio_spec.Step{1}.Boundary.rigid_body{2}.prescribed.ATTR.bc='y';
febio_spec.Step{1}.Boundary.rigid_body{2}.prescribed.ATTR.lc=1;
febio_spec.Step{1}.Boundary.rigid_body{2}.prescribed.VAL=clampCompressiveDisplacement;

febio_spec.Step{1}.Boundary.rigid_body{3}.ATTR.mat=4;
febio_spec.Step{1}.Boundary.rigid_body{3}.fixed{1}.ATTR.bc='x';
febio_spec.Step{1}.Boundary.rigid_body{3}.fixed{2}.ATTR.bc='z';
febio_spec.Step{1}.Boundary.rigid_body{3}.fixed{3}.ATTR.bc='Rx';
febio_spec.Step{1}.Boundary.rigid_body{3}.fixed{4}.ATTR.bc='Ry';
febio_spec.Step{1}.Boundary.rigid_body{3}.fixed{5}.ATTR.bc='Rz';
febio_spec.Step{1}.Boundary.rigid_body{3}.prescribed.ATTR.bc='y';
febio_spec.Step{1}.Boundary.rigid_body{3}.prescribed.ATTR.lc=1;
febio_spec.Step{1}.Boundary.rigid_body{3}.prescribed.VAL=-clampCompressiveDisplacement;

febio_spec.Step{1}.Boundary.rigid_body{4}.ATTR.mat=5;
febio_spec.Step{1}.Boundary.rigid_body{4}.fixed{1}.ATTR.bc='x';
febio_spec.Step{1}.Boundary.rigid_body{4}.fixed{2}.ATTR.bc='z';
febio_spec.Step{1}.Boundary.rigid_body{4}.fixed{3}.ATTR.bc='Rx';
febio_spec.Step{1}.Boundary.rigid_body{4}.fixed{4}.ATTR.bc='Ry';
febio_spec.Step{1}.Boundary.rigid_body{4}.fixed{5}.ATTR.bc='Rz';
febio_spec.Step{1}.Boundary.rigid_body{4}.prescribed.ATTR.bc='y';
febio_spec.Step{1}.Boundary.rigid_body{4}.prescribed.ATTR.lc=1;
febio_spec.Step{1}.Boundary.rigid_body{4}.prescribed.VAL=clampCompressiveDisplacement;

%STEP 2 Tension
%Set 1
febio_spec.Step{2}.Boundary.rigid_body{1}.ATTR.mat=2;
febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='x';
febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='Rx';
febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Ry';
febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Rz';
febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.bc='z';
febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.lc=2;
febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.VAL=clampTensionDisplacement;

febio_spec.Step{2}.Boundary.rigid_body{2}.ATTR.mat=3;
febio_spec.Step{2}.Boundary.rigid_body{2}.fixed{1}.ATTR.bc='x';
febio_spec.Step{2}.Boundary.rigid_body{2}.fixed{2}.ATTR.bc='y';
febio_spec.Step{2}.Boundary.rigid_body{2}.fixed{3}.ATTR.bc='Rx';
febio_spec.Step{2}.Boundary.rigid_body{2}.fixed{4}.ATTR.bc='Ry';
febio_spec.Step{2}.Boundary.rigid_body{2}.fixed{5}.ATTR.bc='Rz';
febio_spec.Step{2}.Boundary.rigid_body{2}.prescribed.ATTR.bc='z';
febio_spec.Step{2}.Boundary.rigid_body{2}.prescribed.ATTR.lc=2;
febio_spec.Step{2}.Boundary.rigid_body{2}.prescribed.VAL=clampTensionDisplacement;

febio_spec.Step{2}.Boundary.rigid_body{3}.ATTR.mat=4;
febio_spec.Step{2}.Boundary.rigid_body{3}.fixed{1}.ATTR.bc='x';
febio_spec.Step{2}.Boundary.rigid_body{3}.fixed{2}.ATTR.bc='y';
febio_spec.Step{2}.Boundary.rigid_body{3}.fixed{3}.ATTR.bc='Rx';
febio_spec.Step{2}.Boundary.rigid_body{3}.fixed{4}.ATTR.bc='Ry';
febio_spec.Step{2}.Boundary.rigid_body{3}.fixed{5}.ATTR.bc='Rz';
febio_spec.Step{2}.Boundary.rigid_body{3}.prescribed.ATTR.bc='z';
febio_spec.Step{2}.Boundary.rigid_body{3}.prescribed.ATTR.lc=2;
febio_spec.Step{2}.Boundary.rigid_body{3}.prescribed.VAL=-clampTensionDisplacement;

febio_spec.Step{2}.Boundary.rigid_body{4}.ATTR.mat=5;
febio_spec.Step{2}.Boundary.rigid_body{4}.fixed{1}.ATTR.bc='x';
febio_spec.Step{2}.Boundary.rigid_body{4}.fixed{2}.ATTR.bc='y';
febio_spec.Step{2}.Boundary.rigid_body{4}.fixed{3}.ATTR.bc='Rx';
febio_spec.Step{2}.Boundary.rigid_body{4}.fixed{4}.ATTR.bc='Ry';
febio_spec.Step{2}.Boundary.rigid_body{4}.fixed{5}.ATTR.bc='Rz';
febio_spec.Step{2}.Boundary.rigid_body{4}.prescribed.ATTR.bc='z';
febio_spec.Step{2}.Boundary.rigid_body{4}.prescribed.ATTR.lc=2;
febio_spec.Step{2}.Boundary.rigid_body{4}.prescribed.VAL=-clampTensionDisplacement;

%Contact section
for q=1:1:4
%     febio_spec.Contact.contact{q}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{q}.ATTR.name;
%     febio_spec.Contact.contact{q}.ATTR.type='sticky';
%     febio_spec.Contact.contact{q}.penalty=100;
%     febio_spec.Contact.contact{q}.laugon=0;
%     febio_spec.Contact.contact{q}.tolerance=0.1;
%     febio_spec.Contact.contact{q}.minaug=0;
%     febio_spec.Contact.contact{q}.maxaug=10;
%     febio_spec.Contact.contact{q}.snap_tol=0;
%     febio_spec.Contact.contact{q}.max_traction=0;
%     febio_spec.Contact.contact{q}.search_tolerance=0.1;
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

%LoadData section
febio_spec.LoadData.loadcurve{1}.ATTR.id=1;
febio_spec.LoadData.loadcurve{1}.ATTR.type='linear';
febio_spec.LoadData.loadcurve{1}.point.VAL=[0 0; 1 1; 2 1];

febio_spec.LoadData.loadcurve{2}.ATTR.id=2;
febio_spec.LoadData.loadcurve{2}.ATTR.type='linear';
febio_spec.LoadData.loadcurve{2}.point.VAL=[0 0; 1 0; 2 1];

%Output section
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_strain;
febio_spec.Output.logfile.element_data{1}.ATTR.data='Ez';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.element_data{1}.VAL=1:1:size(E,1);

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
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]); 
    
    %%
    % Importing element strain energies from a log file
    [~,E_strain,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_strain)); %Element strain energy
    
    %Remove nodal index column
    E_strain=E_strain(:,2:end,:);
    
    %Add initial state i.e. zero energy
    sizImport=size(E_strain); 
    sizImport(3)=sizImport(3)+1;
    E_energy_mat_n=zeros(sizImport);
    E_energy_mat_n(:,:,2:end)=E_strain;
    E_strain=E_energy_mat_n;
    
    [F,CF]=element2patch(E,E_strain(:,:,1));
    CF_V=faceToVertexMeasure(F,V,CF);
    
    %%
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations
    
    axLim=[min(min(V_DEF,[],3),[],1); max(max(V_DEF,[],3),[],1)];
        
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; hold on;%Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1=gpatch(Fb,V_DEF(:,:,end),CF_V,'k',1);
    hp1.FaceColor='interp';
    hp2=gpatch(Fc,V_DEF(:,:,end),'kw','none',0.5);
%     patchNormPlot(Fc,V);
    axisGeom(gca,fontSize);
    colormap(warmcold(250)); colorbar;
    caxis([-tensileStrain tensileStrain]);
    axis(axLim(:)'); %Set axis limits statically

    camlight headlight; axis off;
   
    % Set up animation features
    animStruct.Time=time_mat; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments
        V_def=V+N_disp_mat(:,:,qt); %Current nodal coordinates
        
        [~,CF]=element2patch(E,E_strain(:,:,qt));
        CF_V=faceToVertexMeasure(F,V,CF);
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1 hp1 hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,CF_V,V_def}; %Property values for to set in order to animate
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

