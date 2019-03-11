%% DEMO_febio_0020_vessel_balloon_inflate
% Below is a demonstration for:
%
% * Building geometry for a cylindrical vessel with tetrahedral elements
% * Defining the boundary conditions
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement results

%% Keywords
%
% * febio_spec version 2.5
% * febio, FEBio
% * vessel, cylinder, balloon
% * prescribed displacement
% * contact, sliding
% * hexahedral elements, hex8
% * tube, cylindrical
% * static, solid, multistep
% * hyperelastic, Ogden
% * displacement logfile
% * stress logfile

%%

clear; close all; clc;

%% Plot settings
fontSize=20;
faceAlpha1=0.8;
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
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress
febioLogFileName_strain=[febioFebFileNamePart,'_strain_out.txt']; %Log file name for exporting strain
febioLogFileName_volumeRatio=[febioFebFileNamePart,'_volumeRatio_out.txt']; %Log file name for exporting volume ratio

%Contact parameters
%***
contactInitialOffset=0.1;
contactAlg=2;
switch contactAlg
    case 1
        contactType='sticky';
    case 2
        contactType='facet-to-facet sliding';
    case 3
        contactType='sliding_with_gaps';
        fric_coeff=0;
        fric_penalty=1;
    case 4
        contactType='sliding2';
end

%Specifying geometry parameters vessel (mm)***
pointSpacing=1;
radiusOuter1=6.5/2;
radiusInner1=5.4/2;
radiusOuter2=5.8/2;
radiusInner2=4.7/2;
vesselLength=85;

radiusBalloon=min([radiusInner1 radiusInner2])-contactInitialOffset;
pointSpacingBalloon=pointSpacing/2;
balloonExtraLength=pointSpacing/2;

%Define applied bc's
radialDisplacement=2; %radial displacement after touch
radialDisplacementTotal=radialDisplacement+contactInitialOffset; %Total radial displacement
numSteps=1; %Number of steps

%Material parameter set
c1=1e-3; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
k_factor=500; %Bulk modulus factor
k=c1*k_factor; %Bulk modulus
fomulationType=1; %1=uncoupled, 2=coupled

% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=0; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
runMode='external';

%Contact parameters
contactPenalty=1;

%Visualization parameters
colorLimits_volumeRatio=[0.99 1.01];

%% Creating model boundary polygons
%

nRad=round((2*pi*max([radiusOuter1 radiusOuter2]))/pointSpacing); %Number of radial steps

t=linspace(0,2*pi,nRad)'; %Angles
t=t(1:end-1); %take away last which equals start
v1_Outer=[-(vesselLength/2)*ones(size(t)) radiusOuter1*sin(t) radiusOuter1*cos(t)]; %Circular coordinates

t=linspace(0,2*pi,nRad)'; %Angles
t=t(1:end-1); %take away last which equals start
v2_Outer=[(vesselLength/2)*ones(size(t)) radiusOuter2*sin(t) radiusOuter2*cos(t)]; %Circular coordinates

t=linspace(0,2*pi,nRad)'; %Angles
t=t(1:end-1); %take away last which equals start
v1_Inner=[-(vesselLength/2)*ones(size(t)) radiusInner1*sin(t) radiusInner1*cos(t)]; %Circular coordinates

t=linspace(0,2*pi,nRad)'; %Angles
t=t(1:end-1); %take away last which equals start
v2_Inner=[(vesselLength/2)*ones(size(t)) radiusInner2*sin(t) radiusInner2*cos(t)]; %Circular coordinates

t=linspace(0,2*pi,nRad)'; %Angles
t=t(1:end-1); %take away last which equals start
v1_balloon=[-((vesselLength/2)+balloonExtraLength)*ones(size(t)) radiusBalloon*sin(t) radiusBalloon*cos(t)]; %Circular coordinates
v2_balloon=[ ((vesselLength/2)+balloonExtraLength)*ones(size(t)) radiusBalloon*sin(t) radiusBalloon*cos(t)]; %Circular coordinates

%%
% Plotting model boundary polygons

cFigure;
hold on;
title('Model boundary polygons','FontSize',fontSize);
plotV(v1_Outer,'r.-')
plotV(v1_Inner,'g.-')
plotV(v2_Outer,'b.-')
plotV(v2_Inner,'y.-')
plotV(v1_balloon,'c.-')
plotV(v2_balloon,'k.-')
axisGeom(gca,fontSize);

drawnow;

%% Creating model boundary surfaces

controlStructLoft.numSteps=ceil(vesselLength./pointSpacing);
controlStructLoft.closeLoopOpt=1;
controlStructLoft.patchType='quad';

%Meshing outer surface
[F1,V1,indStart1]=polyLoftLinear(v1_Outer,v2_Outer,controlStructLoft);

%Meshing inner surface
[F2,V2,indStart2]=polyLoftLinear(v1_Inner,v2_Inner,controlStructLoft);

%Compose hexahedral elements
indStart2=indStart2+size(V1,1);
F2=F2+size(V1,1);
E=[F1 F2];
V=[V1;V2];
[FE]=element2patch(E,[],'hex8');

%Meshing balloon surface
[Fs,Vs]=polyLoftLinear(v1_balloon,v2_balloon,controlStructLoft);
Fs=fliplr(Fs); %Invert orientation
[Fs,Vs]=subQuad(Fs,Vs,1);

%%
% Plotting model boundary surfaces

cFigure;
hold on;
title('Model boundary surfaces','FontSize',fontSize);

gpatch(FE,V,'rw','k',0.5);
patchNormPlot(FE,V);

gpatch(F1,V,'bw');
patchNormPlot(F1,V);

gpatch(F2,V,'gw');
patchNormPlot(F2,V);

gpatch(Fs,Vs,'kw');
patchNormPlot(Fs,Vs);

colormap(gca,gjet(4));
icolorbar;

axisGeom(gca,fontSize);
camlight headlight;
drawnow;


%% Joining node sets

Fs=Fs+size(V,1);
V=[V;Vs];

%% Defining the boundary conditions
% The visualization of the model boundary shows colors for each side of the
% cube. These labels can be used to define boundary conditions.

%Define X supported node set
bcSupportList_X=unique([indStart1(:);indStart2(:)]); %Node set part of selected face

% %Define Y supported node set
% bcSupportList_Y=unique(Fb(Cb==4,:)); %Node set part of selected face

%Radial expansion prescribed displacement
bcPrescribeList=(size(V,1)-size(Vs,1)+1):size(V,1);
radialDisplacementStep=radialDisplacementTotal/numSteps; %The radial displacement increment for each step
[th,r,z] = cart2pol(V(bcPrescribeList,2),V(bcPrescribeList,3),V(bcPrescribeList,1));
r2=r+radialDisplacementStep;
V2=V;
[V2(bcPrescribeList,2),V2(bcPrescribeList,3),V2(bcPrescribeList,1)] = pol2cart(th,r2,z);

%%
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists.

hf=cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(FE,V,'kw','k',0.5);
gpatch(Fs,V,'kw','k',0.5);

hl(1)=plotV(V(bcSupportList_X,:),'k.','MarkerSize',markerSize);
hl(2)=gpatch(Fs,V2,'rw','r',0.5);

legend(hl,{'BC X support','BC prescribe 1 step'});

axisGeom(gca,fontSize);
camlight headlight;
drawnow;

%% Define contact

F_contact_master=Fs;
F_contact_slave=F2;

%%
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists.

hf=cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(FE,V,'kw','none',0.5);

h1(1)=gpatch(F_contact_master,V,'gw','g',0.5);
patchNormPlot(F_contact_master,V);
h1(2)=gpatch(F_contact_slave,V,'rw','r',0.5);
patchNormPlot(F_contact_slave,V);

legend(hl,{'Master surface','Slave surface'});

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

%Add template based default settings to proposed control section
[stepStruct.Control]=structComplete(stepStruct.Control,febio_spec.Control,1); %Complement provided with default if missing

%Remove control field (part of template) since step specific control sections are used
febio_spec=rmfield(febio_spec,'Control');

%Material section
if fomulationType==1
    febio_spec.Material.material{1}.ATTR.type='Ogden';
    febio_spec.Material.material{1}.ATTR.id=1;
    febio_spec.Material.material{1}.c1=c1;
    febio_spec.Material.material{1}.m1=m1;
    febio_spec.Material.material{1}.c2=c1;
    febio_spec.Material.material{1}.m2=-m1;
    febio_spec.Material.material{1}.k=k;
    
    febio_spec.Material.material{2}.ATTR.type='Ogden';
    febio_spec.Material.material{2}.ATTR.id=2;
    febio_spec.Material.material{2}.c1=c1;
    febio_spec.Material.material{2}.m1=m1;
    febio_spec.Material.material{2}.c2=c1;
    febio_spec.Material.material{2}.m2=-m1;
    febio_spec.Material.material{2}.k=k;
elseif fomulationType==2
    febio_spec.Material.material{1}.ATTR.type='Ogden unconstrained';
    febio_spec.Material.material{1}.ATTR.id=1;
    febio_spec.Material.material{1}.c1=c1;
    febio_spec.Material.material{1}.m1=m1;
    febio_spec.Material.material{1}.c2=c1;
    febio_spec.Material.material{1}.m2=-m1;
    febio_spec.Material.material{1}.cp=k;
    
    febio_spec.Material.material{2}.ATTR.type='Ogden unconstrained';
    febio_spec.Material.material{2}.ATTR.id=2;
    febio_spec.Material.material{2}.c1=c1;
    febio_spec.Material.material{2}.m1=m1;
    febio_spec.Material.material{2}.c2=c1;
    febio_spec.Material.material{2}.m2=-m1;
    febio_spec.Material.material{2}.cp=k;
end

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='hex8'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set
febio_spec.Geometry.Elements{1}.ATTR.name='Vessel'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=E;

febio_spec.Geometry.Elements{2}.ATTR.type='quad4'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat=2; %material index for this set
febio_spec.Geometry.Elements{2}.ATTR.name='Balloon'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id=size(E,1)+(1:1:size(Fs,1))'; %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL=Fs;

% -> NodeSets
febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList_X';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcSupportList_X(:);

febio_spec.Geometry.NodeSet{2}.ATTR.name='bcPrescribeList';
febio_spec.Geometry.NodeSet{2}.node.ATTR.id=bcPrescribeList(:);

% -> Surfaces
febio_spec.Geometry.Surface{1}.ATTR.name='contact_master';
febio_spec.Geometry.Surface{1}.quad4.ATTR.lid=(1:1:size(F_contact_master,1))';
febio_spec.Geometry.Surface{1}.quad4.VAL=F_contact_master;

febio_spec.Geometry.Surface{2}.ATTR.name='contact_slave';
febio_spec.Geometry.Surface{2}.quad4.ATTR.lid=(1:1:size(F_contact_slave,1))';
febio_spec.Geometry.Surface{2}.quad4.VAL=F_contact_slave;

% -> Surface pairs
febio_spec.Geometry.SurfacePair{1}.ATTR.name='Contact1';
febio_spec.Geometry.SurfacePair{1}.master.ATTR.surface=febio_spec.Geometry.Surface{1}.ATTR.name;
febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;

%MeshData section
febio_spec.MeshData.ElementData{1}.ATTR.var='shell thickness';
febio_spec.MeshData.ElementData{1}.ATTR.elem_set=febio_spec.Geometry.Elements{2}.ATTR.name;
febio_spec.MeshData.ElementData{1}.elem.ATTR.lid=(1:size(Fs,1))';
febio_spec.MeshData.ElementData{1}.elem.VAL=0.1*ones(size(Fs));

%Create steps
[th,r2,z] = cart2pol(V(bcPrescribeList,2),V(bcPrescribeList,3),V(bcPrescribeList,1));
nodeSetName=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.MeshData.NodeData=[];%Initialize so we can use end+1 indexing
bcNames={'x','y','z'};

if numSteps==1
    
    febio_spec.Control.analysis.ATTR.type='static';
    febio_spec.Control.time_steps=numTimeSteps;
    febio_spec.Control.step_size=1/numTimeSteps;
    febio_spec.Control.time_stepper.dtmin=dtmin;
    febio_spec.Control.time_stepper.dtmax=dtmax;
    febio_spec.Control.time_stepper.max_retries=max_retries;
    febio_spec.Control.time_stepper.opt_iter=opt_iter;
    febio_spec.Control.max_refs=max_refs;
    febio_spec.Control.max_ups=max_ups;
    
    %Define prescribed displacements
    bcPrescribeMagnitudesStep=V2(bcPrescribeList,:)-V(bcPrescribeList,:);
    
    %Define mesh data and prescribed displacements
    for q_dir=1:1:3 %Loop over coordinates dimensions
        
        %Define mesh data for displacement increments
        c=numel(febio_spec.MeshData.NodeData)+1; %Current step index
        febio_spec.MeshData.NodeData{c}.ATTR.name=['displacement_',bcNames{q_dir},'_1'];
        febio_spec.MeshData.NodeData{c}.ATTR.node_set=nodeSetName;
        febio_spec.MeshData.NodeData{c}.node.ATTR.lid=(1:1:numel(bcPrescribeList))';
        febio_spec.MeshData.NodeData{c}.node.VAL=bcPrescribeMagnitudesStep(:,q_dir);
        
        %Define prescribed displacements
        febio_spec.Boundary.prescribe{q_dir}.ATTR.bc=bcNames{q_dir};
        febio_spec.Boundary.prescribe{q_dir}.ATTR.relative=1;
        febio_spec.Boundary.prescribe{q_dir}.ATTR.node_set=nodeSetName;
        febio_spec.Boundary.prescribe{q_dir}.scale.ATTR.lc=1;
        febio_spec.Boundary.prescribe{q_dir}.scale.VAL=1;
        febio_spec.Boundary.prescribe{q_dir}.relative=1;
        febio_spec.Boundary.prescribe{q_dir}.value.ATTR.node_data=febio_spec.MeshData.NodeData{c}.ATTR.name;
    end
else
    V2n=V;
    V2=V;
    for q=1:1:numSteps
        %Step specific control section
        febio_spec.Step{q}.ATTR.id=q;
        febio_spec.Step{q}.Control=stepStruct.Control;
        
        %Offset coordinates
        r2=r2+radialDisplacementStep;
        V2n=V2;
        [V2n(bcPrescribeList,2),V2n(bcPrescribeList,3),V2n(bcPrescribeList,1)] = pol2cart(th,r2,z); %The current set
        
        %Define prescribed displacements
        bcPrescribeMagnitudesStep=V2n(bcPrescribeList,:)-V2(bcPrescribeList,:);
        V2=V2n;
        
        %Define mesh data and prescribed displacements
        for q_dir=1:1:3 %Loop over coordinates dimensions
            
            %Define mesh data for displacement increments
            c=numel(febio_spec.MeshData.NodeData)+1; %Current step index
            febio_spec.MeshData.NodeData{c}.ATTR.name=['displacement_',bcNames{q_dir},'_step_',num2str(q)];
            febio_spec.MeshData.NodeData{c}.ATTR.node_set=nodeSetName;
            febio_spec.MeshData.NodeData{c}.node.ATTR.lid=(1:1:numel(bcPrescribeList))';
            febio_spec.MeshData.NodeData{c}.node.VAL=bcPrescribeMagnitudesStep(:,q_dir);
            
            %Define prescribed displacements
            febio_spec.Step{q}.Boundary.prescribe{q_dir}.ATTR.bc=bcNames{q_dir};
            febio_spec.Step{q}.Boundary.prescribe{q_dir}.ATTR.relative=1;
            febio_spec.Step{q}.Boundary.prescribe{q_dir}.ATTR.node_set=nodeSetName;
            febio_spec.Step{q}.Boundary.prescribe{q_dir}.scale.ATTR.lc=1;
            febio_spec.Step{q}.Boundary.prescribe{q_dir}.scale.VAL=1;
            febio_spec.Step{q}.Boundary.prescribe{q_dir}.relative=1;
            febio_spec.Step{q}.Boundary.prescribe{q_dir}.value.ATTR.node_data=febio_spec.MeshData.NodeData{c}.ATTR.name;
        end
        
    end
    
end
%Boundary condition section
% -> Fix boundary conditions
febio_spec.Boundary.fix{1}.ATTR.bc='x';
febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
% febio_spec.Boundary.fix{2}.ATTR.bc='y';
% febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
% febio_spec.Boundary.fix{3}.ATTR.bc='z';
% febio_spec.Boundary.fix{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;

%Contact section
switch contactType
    case 'sticky'
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='sticky';
        febio_spec.Contact.contact{1}.penalty=contactPenalty;
        febio_spec.Contact.contact{1}.laugon=0;
        febio_spec.Contact.contact{1}.tolerance=0.1;
        febio_spec.Contact.contact{1}.minaug=0;
        febio_spec.Contact.contact{1}.maxaug=10;
        febio_spec.Contact.contact{1}.snap_tol=0;
        febio_spec.Contact.contact{1}.max_traction=0;
        febio_spec.Contact.contact{1}.search_tolerance=0.1;
    case 'facet-to-facet sliding'
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='facet-to-facet sliding';
        febio_spec.Contact.contact{1}.penalty=contactPenalty;
        febio_spec.Contact.contact{1}.auto_penalty=1;
        febio_spec.Contact.contact{1}.two_pass=0;
        febio_spec.Contact.contact{1}.laugon=0;
        febio_spec.Contact.contact{1}.tolerance=0.1;
        febio_spec.Contact.contact{1}.gaptol=0;
        febio_spec.Contact.contact{1}.minaug=0;
        febio_spec.Contact.contact{1}.maxaug=10;
        febio_spec.Contact.contact{1}.search_tol=0.01;
        febio_spec.Contact.contact{1}.search_radius=pointSpacing/10;
    case 'sliding_with_gaps'
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='sliding_with_gaps';
        febio_spec.Contact.contact{1}.penalty=contactPenalty;
        febio_spec.Contact.contact{1}.auto_penalty=1;
        febio_spec.Contact.contact{1}.two_pass=0;
        febio_spec.Contact.contact{1}.laugon=0;
        febio_spec.Contact.contact{1}.tolerance=0.1;
        febio_spec.Contact.contact{1}.gaptol=0;
        febio_spec.Contact.contact{1}.minaug=0;
        febio_spec.Contact.contact{1}.maxaug=10;
        febio_spec.Contact.contact{1}.fric_coeff=fric_coeff;
        febio_spec.Contact.contact{1}.fric_penalty=fric_penalty;
        febio_spec.Contact.contact{1}.ktmult=1;
        febio_spec.Contact.contact{1}.seg_up=0;
        febio_spec.Contact.contact{1}.search_tol=0.01;
    case 'sliding2'
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='sliding2';
        febio_spec.Contact.contact{1}.penalty=contactPenalty;
        febio_spec.Contact.contact{1}.auto_penalty=1;
        febio_spec.Contact.contact{1}.two_pass=0;
        febio_spec.Contact.contact{1}.laugon=0;
        febio_spec.Contact.contact{1}.tolerance=0.1;
        febio_spec.Contact.contact{1}.gaptol=0;
        febio_spec.Contact.contact{1}.symmetric_stiffness=0;
        febio_spec.Contact.contact{1}.search_tol=0.01;
        febio_spec.Contact.contact{1}.search_radius=pointSpacing/2;
end

% LoadData section
febio_spec.LoadData.loadcurve{1}.ATTR.id=1;
febio_spec.LoadData.loadcurve{1}.ATTR.type='linear';
febio_spec.LoadData.loadcurve{1}.point.VAL=[0 0; 1 1;];

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
febio_spec.Output.logfile.element_data{1}.VAL=1:size(E,1);

febio_spec.Output.logfile.element_data{2}.ATTR.file=febioLogFileName_strain;
febio_spec.Output.logfile.element_data{2}.ATTR.data='E1;E2;E3';
febio_spec.Output.logfile.element_data{2}.ATTR.delim=',';
febio_spec.Output.logfile.element_data{2}.VAL=1:size(E,1);

febio_spec.Output.logfile.element_data{3}.ATTR.file=febioLogFileName_volumeRatio;
febio_spec.Output.logfile.element_data{3}.ATTR.data='J';
febio_spec.Output.logfile.element_data{3}.ATTR.delim=',';
febio_spec.Output.logfile.element_data{3}.VAL=1:size(E,1);

%% Quick viewing of the FEBio input file structure
% The |febView| function can be used to view the xml structure in a MATLAB
% figure window.

%%
% |febView(febio_spec); %Viewing the febio file|
% febView(febio_spec)

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
febioAnalysis.maxLogCheckTime=30; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results

if runFlag==1 | runFlag==0 %i.e. a succesful run
    
    % Importing nodal displacements from a log file
    [time_mat, N_disp_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp)); %Nodal displacements
    time_mat=[0; time_mat(:)]; %Time
    
    N_disp_mat=N_disp_mat(:,2:end,:);
    sizImport=size(N_disp_mat);
    sizImport(3)=sizImport(3)+1;
    N_disp_mat_n=zeros(sizImport);
    N_disp_mat_n(:,:,2:end)=N_disp_mat;
    N_disp_mat=N_disp_mat_n;
    
    DN_MAG=sqrt(sum(N_disp_mat.^2,2));
    DN=N_disp_mat(:,:,end);
    DN_magnitude=sqrt(sum(DN(:,3).^2,2));
    V_def=V+DN;
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
    X_DEF=V_DEF(:,1,:);
    Y_DEF=V_DEF(:,2,:);
    Z_DEF=V_DEF(:,3,:);
    %     [CF]=vertexToFaceMeasure(Fb,DN_magnitude);
    
    for q=1:1:numel(febio_spec.Output.logfile.element_data)
        
        elementDataFileName=febio_spec.Output.logfile.element_data{q}.ATTR.file;
        
        % Importing element data from log file
        [~,E_data,~]=importFEBio_logfile(fullfile(savePath,elementDataFileName)); %Element data
        
        %Remove nodal index column
        E_data=E_data(:,2:end,:);
        
        %Add initial state
        sizImport=size(E_data);
        sizImport(3)=sizImport(3)+1;
        if strcmp('J',febio_spec.Output.logfile.element_data{q}.ATTR.data)
            E_data_mat_n=ones(sizImport);
        else
            E_data_mat_n=zeros(sizImport);
        end
        E_data_mat_n(:,:,2:end)=E_data;
        E_data=E_data_mat_n;
        
        for qd=1:1:size(E_data,2)
            E_data_now=E_data(:,qd,:);
            
            [F,CF]=element2patch(E,E_data_now(:,:,1));
            
            % Plotting the simulated results using |anim8| to visualize and animate
            % deformations
            
            % Create basic view and store graphics handle to initiate animation
            hf=cFigure; %Open figure
            gtitle([febioFebFileNamePart,', data: ',febio_spec.Output.logfile.element_data{q}.ATTR.data,', index: ',num2str(qd)]);
            
            hold on;
            hp1=gpatch(F,V_def,CF,'k',1); %Add graphics object to animate
            hp2=gpatch(Fs,V_def,'kw','k',1); %Add graphics object to animate
            %     gpatch(FE,V,0.5*ones(1,3),'k',0.25); %A static graphics object
            
            colormap(gca,gjet(250)); hc=colorbar;
            
            caxis([min(E_data_now(:)) max(E_data_now(:))]);
            axisGeom(gca,fontSize);
            axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);
            axis manual;
            camlight headlight;
            drawnow;
            
            % Set up animation features
            animStruct.Time=time_mat; %The time vector
            for qt=1:1:size(N_disp_mat,3) %Loop over time increments
                DN=N_disp_mat(:,:,qt); %Current displacement
                DN_magnitude=sqrt(sum(DN.^2,2)); %Current displacement magnitude
                V_def=V+DN; %Current nodal coordinates
                [~,CF]=element2patch(E,E_data_now(:,:,qt));
                
                %Set entries in animation structure
                animStruct.Handles{qt}=[hp1 hp1 hp2]; %Handles of objects to animate
                animStruct.Props{qt}={'Vertices','CData','Vertices',}; %Properties of objects to animate
                animStruct.Set{qt}={V_def,CF,V_def}; %Property values for to set in order to animate
            end
            anim8(hf,animStruct); %Initiate animation feature
            drawnow;
            
        end
    end
    
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
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
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
