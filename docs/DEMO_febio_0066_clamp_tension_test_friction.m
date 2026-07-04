%% DEMO_febio_0066_clamp_tension_test_friction
% Below is a demonstration for: 
% 1) The creation of an FEBio model for clamped tensile testing
% 2) The use of multiple steps
% 4) Running an FEBio job with MATLAB
% 5) Importing FEBio results into MATLAB

%% Keywords
%
% * febio_spec version 4.0
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
pointSpacing=1; 
sampleWidth=10;
sampleThickness=2; 
sampleClampedHeight=sampleWidth;
sampleGripGripHeight=sampleWidth.*2;
appliedLinearStrain=0.3;
clampCompressiveLinearStrain=0.3;

pointSpacingClamp=pointSpacing; 
clampExtendAmount=sampleWidth/2;

initialSpacing=0;%pointSpacing/10;

%Initial material parameter set
c1=1e-3;
m1=2;
k_factor=100;
k=c1*k_factor;

% FEA control settings
numTimeSteps=20; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=25; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
min_residual=1e-20;
symmetric_stiffness=0;

runMode='external'; %'internal' or 'external'

%Contact parameters
contactPenalty=10;
contactType='sliding-elastic'; % 'sticky'
laugon=0;
minaug=1;
maxaug=10;
fric_coeff=1;

%% Computing derived parameters
sampleHeight=sampleGripGripHeight+2*sampleClampedHeight;

numElementsWidthStrip=ceil(sampleWidth/pointSpacing);
numElementsWidthStrip=numElementsWidthStrip+iseven(numElementsWidthStrip); %Force uneven so there is a middle element
numElementsThicknessStrip=ceil(round(sampleThickness/pointSpacing));
numElementsHeightStrip=ceil(round(sampleHeight/pointSpacing));
numElementsHeightStrip=numElementsHeightStrip+iseven(numElementsHeightStrip); %Force uneven so there is a middle element

clampWidth=sampleWidth+(2*clampExtendAmount);
clampHeight=sampleClampedHeight+(clampExtendAmount);
numElementsWidthClamp=ceil(clampWidth/pointSpacingClamp);
numElementsHeightClamp=ceil(clampHeight/pointSpacingClamp);

clampCompressiveDisplacement=((sampleThickness.*clampCompressiveLinearStrain)/2)+initialSpacing;
clampTensionDisplacement=(sampleGripGripHeight.*appliedLinearStrain);

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
patchNormPlot(Fb1,V);

gpatch(Fb2,V,'kw','k',1);
patchNormPlot(Fb2,V);

gpatch(Fc1,V,'rw','k',1);
patchNormPlot(Fc1,V);

gpatch(Fc2,V,'gw','k',1);
patchNormPlot(Fc2,V);

gpatch(Fc3,V,'bw','k',1);
patchNormPlot(Fc3,V);

gpatch(Fc4,V,'yw','k',1);
patchNormPlot(Fc4,V);

axisGeom(gca,fontSize);
camlight headlight;
drawnow; 

%% Get logic for middle elements (to help investigate strain in the middle)

searchTol=pointSpacing/1000; 
VE=patchCentre(E,V);
logicMiddleElements= VE(:,3)<searchTol & VE(:,3)>-searchTol;

[FM]=element2patch(E(logicMiddleElements,:),[],'hex8');

%%
cFigure; hold on;
gpatch(Fb,V,'kw','k',0.25);

gpatch(FM,V,'rw','k',1);

% hp(1)=plotV(V(bcPrescribeList1,:),'r.','MarkerSize',25);
% legend(hp,{'Node set 1','Node set 2','Node set 3','Node set 4'})
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
stepStruct.Control.solver.symmetric_stiffness=symmetric_stiffness;

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
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.ATTR.type='Ogden';
febio_spec.Material.material{1}.c1=c1;
febio_spec.Material.material{1}.m1=m1;
febio_spec.Material.material{1}.c2=c1;
febio_spec.Material.material{1}.m2=-m1;
febio_spec.Material.material{1}.k=k;

materialName2='Material2';
febio_spec.Material.material{2}.ATTR.name=materialName2;
febio_spec.Material.material{2}.ATTR.type='rigid body';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.density=1;
febio_spec.Material.material{2}.center_of_mass=mean(Vc1,1);

materialName3='Material3';
febio_spec.Material.material{3}.ATTR.name=materialName3;
febio_spec.Material.material{3}.ATTR.type='rigid body';
febio_spec.Material.material{3}.ATTR.id=3;
febio_spec.Material.material{3}.density=1;
febio_spec.Material.material{3}.center_of_mass=mean(Vc2,1);

materialName4='Material4';
febio_spec.Material.material{4}.ATTR.name=materialName4;
febio_spec.Material.material{4}.ATTR.type='rigid body';
febio_spec.Material.material{4}.ATTR.id=4;
febio_spec.Material.material{4}.density=1;
febio_spec.Material.material{4}.center_of_mass=mean(Vc3,1);

materialName5='Material5';
febio_spec.Material.material{5}.ATTR.name=materialName5;
febio_spec.Material.material{5}.ATTR.type='rigid body';
febio_spec.Material.material{5}.ATTR.id=5;
febio_spec.Material.material{5}.density=1;
febio_spec.Material.material{5}.center_of_mass=mean(Vc4,1);

%Geometry section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=mcol(1:size(V,1)); %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1_strip';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of the element set
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type of this set
febio_spec.Mesh.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Mesh.Elements{1}.elem.VAL=E;
febio_spec.Mesh.Elements{1}.elem.ATTR.id=mcol(1:size(E,1)); %Element material id's

partName2='Part2_Grip1';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of the element set
febio_spec.Mesh.Elements{2}.ATTR.type='quad4'; %Element type of this set
febio_spec.Mesh.Elements{2}.ATTR.mat=2; %material index for this set 
febio_spec.Mesh.Elements{2}.elem.VAL=Fc1;
febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E,1)+mcol(1:size(Fc1,1)); %Element material id's

partName3='Part3_Grip2';
febio_spec.Mesh.Elements{3}.ATTR.name=partName3; %Name of the element set
febio_spec.Mesh.Elements{3}.ATTR.type='quad4'; %Element type of this set
febio_spec.Mesh.Elements{3}.ATTR.mat=3; %material index for this set 
febio_spec.Mesh.Elements{3}.elem.VAL=Fc2;
febio_spec.Mesh.Elements{3}.elem.ATTR.id=size(E,1)+size(Fc1,1)+mcol(1:size(Fc2,1)); %Element material id's

partName4='Part4_Grip3';
febio_spec.Mesh.Elements{4}.ATTR.name=partName4; %Name of the element set
febio_spec.Mesh.Elements{4}.ATTR.type='quad4'; %Element type of this set
febio_spec.Mesh.Elements{4}.ATTR.mat=4; %material index for this set 
febio_spec.Mesh.Elements{4}.elem.VAL=Fc3;
febio_spec.Mesh.Elements{4}.elem.ATTR.id=size(E,1)+2*size(Fc1,1)+mcol(1:size(Fc3,1)); %Element material id's

partName5='Part5_Grip4';
febio_spec.Mesh.Elements{5}.ATTR.name=partName5; %Name of the element set
febio_spec.Mesh.Elements{5}.ATTR.type='quad4'; %Element type of this set
febio_spec.Mesh.Elements{5}.ATTR.mat=5; %material index for this set 
febio_spec.Mesh.Elements{5}.elem.VAL=Fc4;
febio_spec.Mesh.Elements{5}.elem.ATTR.id=size(E,1)+3*size(Fc1,1)+mcol(1:size(Fc4,1)); %Element material id's

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

% -> Surfaces
febio_spec.Mesh.Surface{1}.ATTR.name='contact_secondary1';
febio_spec.Mesh.Surface{1}.quad4.ATTR.id=mcol(1:1:size(Fc1,1));
febio_spec.Mesh.Surface{1}.quad4.VAL=Fc1;

febio_spec.Mesh.Surface{2}.ATTR.name='contact_secondary2';
febio_spec.Mesh.Surface{2}.quad4.ATTR.id=mcol(1:1:size(Fc2,1));
febio_spec.Mesh.Surface{2}.quad4.VAL=Fc2;

febio_spec.Mesh.Surface{3}.ATTR.name='contact_secondary3';
febio_spec.Mesh.Surface{3}.quad4.ATTR.id=mcol(1:1:size(Fc3,1));
febio_spec.Mesh.Surface{3}.quad4.VAL=Fc3;

febio_spec.Mesh.Surface{4}.ATTR.name='contact_secondary4';
febio_spec.Mesh.Surface{4}.quad4.ATTR.id=mcol(1:1:size(Fc4,1));
febio_spec.Mesh.Surface{4}.quad4.VAL=Fc4;

febio_spec.Mesh.Surface{5}.ATTR.name='contact_primary1';
febio_spec.Mesh.Surface{5}.quad4.ATTR.id=mcol(1:1:size(Fb1,1));
febio_spec.Mesh.Surface{5}.quad4.VAL=Fb1;

febio_spec.Mesh.Surface{6}.ATTR.name='contact_primary2';
febio_spec.Mesh.Surface{6}.quad4.ATTR.id=mcol(1:1:size(Fb2,1));
febio_spec.Mesh.Surface{6}.quad4.VAL=Fb2;

% -> Surface pairs
for q=1:1:4
    if iseven(q)
        indStrip=5;
    else
        indStrip=6;
    end    
    febio_spec.Mesh.SurfacePair{q}.ATTR.name=['Contact_',num2str(q)];
    febio_spec.Mesh.SurfacePair{q}.primary=febio_spec.Mesh.Surface{indStrip}.ATTR.name;
    febio_spec.Mesh.SurfacePair{q}.secondary=febio_spec.Mesh.Surface{q}.ATTR.name;        
end

%Rigid section 

rbIds=[2 3 4 5]; %Rigid body ids
s=[-1 1 -1 1]; %sign/multiplier
bcId=0;
for q=1:1:numel(rbIds)
    
    bcId=bcId+1;

    % ->Rigid body fix boundary conditions
    febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.ATTR.name=['RigidFix_rb',num2str(rbIds(q)),'_01'];
    febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.ATTR.type='rigid_fixed';
    febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.rb=rbIds(q);
    febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.Rx_dof=1;
    %febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.Ry_dof=1;
    febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.Rz_dof=1;
    febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.Ru_dof=1;
    febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.Rv_dof=1;
    febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.Rw_dof=1;

    bcId=bcId+1;

    % ->Rigid body prescribe boundary conditions
    febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.ATTR.name=['RigidPrescribe_rb',num2str(rbIds(q)),'_01'];
    febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.ATTR.type='rigid_displacement';
    febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.rb=rbIds(q);
    febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.dof='y';
    febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.value.ATTR.lc=1;
    febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.value.VAL=s(q).*clampCompressiveDisplacement;
    febio_spec.Step.step{1}.Rigid.rigid_bc{bcId}.relative=0;

end

rbIds=[2 3 4 5]; %Rigid body ids
s=[1 1 0 0]; %sign/multiplier
bcId=0;
for q=1:1:numel(rbIds)
    
    bcId=bcId+1;

    % ->Rigid body fix boundary conditions
    febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.ATTR.name=['RigidFix_rb',num2str(rbIds(q)),'_02'];
    febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.ATTR.type='rigid_fixed';
    febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.rb=rbIds(q);
    febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.Rx_dof=1;
    febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.Ry_dof=1;
    %febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.Rz_dof=1;
    febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.Ru_dof=1;
    febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.Rv_dof=1;
    febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.Rw_dof=1;

    bcId=bcId+1;

    % ->Rigid body prescribe boundary conditions
    febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.ATTR.name=['RigidPrescribe_rb',num2str(rbIds(q)),'_02'];
    febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.ATTR.type='rigid_displacement';
    febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.rb=rbIds(q);
    febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.dof='z';
    febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.value.ATTR.lc=2;
    febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.value.VAL=s(q).*clampTensionDisplacement;
    febio_spec.Step.step{2}.Rigid.rigid_bc{bcId}.relative=1;

end

%Contact section
for q=1:1:4
    %Contact section
    switch contactType
        case 'sliding-elastic'
            febio_spec.Contact.contact{q}.ATTR.type=contactType;
            febio_spec.Contact.contact{q}.ATTR.surface_pair=febio_spec.Mesh.SurfacePair{q}.ATTR.name;
            febio_spec.Contact.contact{q}.two_pass=0;
            febio_spec.Contact.contact{q}.laugon=laugon;
            febio_spec.Contact.contact{q}.tolerance=0.2;
            febio_spec.Contact.contact{q}.gaptol=0;
            febio_spec.Contact.contact{q}.minaug=minaug;
            febio_spec.Contact.contact{q}.maxaug=maxaug;
            febio_spec.Contact.contact{q}.search_tol=0.01;
            febio_spec.Contact.contact{q}.search_radius=0.1*sqrt(sum((max(V,[],1)-min(V,[],1)).^2,2));
            febio_spec.Contact.contact{q}.symmetric_stiffness=0;
            febio_spec.Contact.contact{q}.auto_penalty=1;
            febio_spec.Contact.contact{q}.update_penalty=1;
            febio_spec.Contact.contact{q}.penalty=contactPenalty;
            febio_spec.Contact.contact{q}.fric_coeff=fric_coeff;
        case 'sticky'
            febio_spec.Contact.contact{q}.ATTR.type=contactType;
            febio_spec.Contact.contact{q}.ATTR.surface_pair=febio_spec.Mesh.SurfacePair{q}.ATTR.name;
                        
            febio_spec.Contact.contact{q}.penalty=contactPenalty;              

            febio_spec.Contact.contact{q}.laugon=laugon;
            febio_spec.Contact.contact{q}.minaug=minaug;
            febio_spec.Contact.contact{q}.maxaug=maxaug;
            febio_spec.Contact.contact{q}.tolerance=0.2;

            febio_spec.Contact.contact{q}.search_tolerance=0.01;
            febio_spec.Contact.contact{q}.max_traction=0;
            febio_spec.Contact.contact{q}.snap_tol=0.01;   
    end
end

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
febio_spec.Output.logfile.element_data{1}.VAL=1:1:size(E,1);

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
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; hold on;%Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('$E_{zz}$','Interpreter','Latex')
    hp1=gpatch(Fb,V_DEF(:,:,end),CV,'k',1);
    hp1.FaceColor='interp';
    hp2=gpatch(Fc,V_DEF(:,:,end),'w','none',0.5);
    
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
        animStruct.Handles{qt}=[hp1 hp1 hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),CV,V_DEF(:,:,qt)}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct); %Initiate animation feature
    gdrawnow;

    %%
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; 
    subplot(1,2,1); hold on;
    title('$E_{zz}$','Interpreter','Latex');
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1=gpatch(Fb,V_DEF(:,:,end),CV,'none',1,0.5); hp1.FaceColor='interp';
    hp2=gpatch(Fc,V_DEF(:,:,end),'w','none',0.5);
    hp3=gpatch(FM,V_DEF(:,:,end),CV,'k',1,2); hp3.FaceColor='interp';
    
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
        % [F,CF]=element2patch(E,E_strain(:,:,qt));
        % CV=faceToVertexMeasure(F,V,CF);
        CV=faceToVertexMeasure(E,V,E_strain(:,:,qt));

        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1 hp1 hp2 hp3 hp3 hp4 hp4 hp5 hp5]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices','Vertices','CData','XData','YData','XData','YData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),CV,V_DEF(:,:,qt),V_DEF(:,:,qt),CV,timeVec(qt),E_strain_middle_mean(qt),timeVec(qt),gripStrain(qt)}; %Property values for to set in order to animate
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
