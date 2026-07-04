%% DEMO_febio_0072_pneunet_actuator_set_disp_01
% Below is a demonstration for:
% 
% * Building geometry for a simple pneunet actuator
% * Defining the boundary conditions 
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement and stress results

%% Keywords
%
% * febio_spec version 4.0
% * febio, FEBio
% * pressure loading
% * hexahedral elements, hex8
% * pneunet actuator
% * soft robotic
% * static, solid
% * hyperelastic, Ogden
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

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

%Load
appliedPressure1=0.15; 
appliedPressure2=appliedPressure1/30; 

%Define applied displacement perturbation
prescribedDisplacement_X=2;

%Material parameter set
c1=1; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
k_factor=100; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus

c2=c1*2; %Shear-modulus-like parameter
m2=2; %Material parameter setting degree of non-linearity
k2=c2*k_factor; %Bulk modulus

% FEA control settings
numTimeSteps=25; %Number of time steps desired
opt_iter=25; %Optimum number of iterations
max_refs=opt_iter*2; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=(1/numTimeSteps)*4; %Maximum time step size

runMode='external';%'internal';

%%

actuatorLengthDesired=100;

pointSpacing=2;
periodSize=10; 
numPeriods=ceil(actuatorLengthDesired/periodSize); 

heightFactor=1; %To scale height wrt period size
widthFactor=1; %To scale width wrt period size

numElementsPeriod=ceil(periodSize./pointSpacing);
if numElementsPeriod<4
    numElementsPeriod=4; 
end

numElementsLength=((numPeriods-1)*numElementsPeriod)+(numElementsPeriod-1);
modelLength=numElementsLength.*pointSpacing;
modelHeight_X=periodSize*heightFactor;
modelWidth_Y=periodSize*widthFactor;

numElementsHeight_X=ceil(modelHeight_X./pointSpacing);
if numElementsHeight_X<6
    numElementsHeight_X=6; 
end

numElementsWidth_Y=ceil(modelWidth_Y./pointSpacing);
if numElementsWidth_Y<1 
    numElementsWidth_Y=1; 
end

boxDim=[modelHeight_X modelWidth_Y modelLength];
boxEl=[numElementsHeight_X numElementsWidth_Y numElementsLength];

[meshStruct]=hexMeshBox(boxDim,boxEl);

E_bar=meshStruct.E;
V_bar=meshStruct.V;
F_bar=meshStruct.F;
Fb_bar=meshStruct.Fb;
Cb_bar=meshStruct.faceBoundaryMarker;

VE_bar=patchCentre(E_bar,V_bar); 

CZ=VE_bar(:,3);
CZ=CZ-min(CZ); 
CZ=CZ./max(CZ); 
CZ=round((CZ.*(numElementsLength-1)))+1;

CW=VE_bar(:,1);
CW=CW-min(CW); 
CW=CW./max(CW); 
CW=round((CW.*(numElementsHeight_X-1)))+1;

CD=rem(CZ,numElementsPeriod); 

logicKeep1=~(CD==0 & CW>3);

E1=E_bar(logicKeep1,:);
F1=element2patch(E1);
[indBoundary1]=tesBoundary(F1,V_bar);

logicKeep2=any(ismember(E1,F1(indBoundary1,:)),2);
F2=element2patch(E1(logicKeep2,:));
[indBoundary2]=tesBoundary(F2,V_bar);

Fb=F2(indBoundary2,:);
Cb=7*ones(size(Fb,1),1);

for q=1:1:6
    F_Cb1=Fb_bar(Cb_bar==q,:);
    logicNow=all(ismember(Fb,F_Cb1),2);
    Cb(logicNow)=q;
end
Cb(~any(ismember(Fb,F1(indBoundary1,:)),2))=0;

%% Removed unused nodes and clean up index matrices
 
[E,V,indFix2]=patchCleanUnused(E1(logicKeep2,:),V_bar);
V(:,1)=V(:,1)-min(V(:,1));
Fb=indFix2(Fb);
F=indFix2(F2);

V2=V;
V2(:,1)=-V2(:,1);
V=[V;V2];

Fb=[Fb;fliplr(Fb)+size(V2,1)];
F=[F;fliplr(F)+size(V2,1)];
cShift=max(Cb(:))+1; 
Cb=[Cb;Cb+cShift];

indicesInnerSurface=unique(Fb(Cb==1,:));
logicElementsInner=any(ismember(E,indicesInnerSurface),2);

E1=E;
E2=E(:,[5 6 7 8 1 2 3 4]);
E2=E2(~logicElementsInner,:);

E=[E1;E2+size(V2,1)];
C=[ones(size(E1,1),1); 2*ones(size(E2,1),1);];
[F,V,~,ind2]=mergeVertices(F,V);
Fb=ind2(Fb);
E=ind2(E);
indicesInnerSurface=ind2(indicesInnerSurface);

actuatorLength=max(V(:,3))-min(V(:,3));

%% Defining the boundary conditions
% The visualization of the model boundary shows colors for each side of the
% disc. These labels can be used to define boundary conditions. 

%Define supported node sets
bcSupportList=unique(Fb(ismember(Cb,[5 5+cShift]),:)); %Node set part of selected face

bcPrescribeList=unique(Fb(ismember(Cb,[6 6+cShift]),:)); %Node set part of selected face

%Get pressure faces
F_pressure1=Fb(Cb==0,:); 
F_pressure2=Fb(Cb==cShift,:); 

%%


ind2=unique(E(C==2,:));
d=max(V(ind2,1));
V(ind2,1)=V(ind2,1)-d;
V(:,1)=V(:,1)+d/2;

% [E,V,indFix]=patchCleanUnused(E,V);
% bcSupportList=indFix(bcSupportList);
% bcPrescribeList=indFix(bcPrescribeList);
% F_pressure1=indFix(F_pressure1);
% F_pressure2=indFix(F_pressure2);

[F,CF]=element2patch(E,C,'hex8');

[F,V,~,indFix]=mergeVertices(F,V);
E=indFix(E);
bcSupportList=indFix(bcSupportList);
bcPrescribeList=indFix(bcPrescribeList);
F_pressure1=indFix(F_pressure1);
F_pressure2=indFix(F_pressure2);

indb=tesBoundary(F,V);
Fb=F(indb,:);

%%

cFigure; hold on;
gpatch(F,V,CF,'k',0.5);
% plotV(V(indicesInnerSurface,:),'k.','MarkerSize',markerSize);
axisGeom; 
colormap(turbo(250)); icolorbar; 
camlight headlight; 
gdrawnow; 

%% Find top surface elements

logicTopSurface=all(ismember(Fb,bcPrescribeList),2);
F_top=Fb(logicTopSurface,:);
center_of_mass=mean(V(unique(F_top(:)),:),1);

%% 
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists. 

hf=cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,V,'w','none',0.5);

hl(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',markerSize);
hl(2)=gpatch(F_pressure1,V,'r','k',1);
hl(3)=gpatch(F_pressure2,V,'b','k',1);
hl(4)=gpatch(F_top,V,'g','k',1);

patchNormPlot(F_pressure1,V);
legend(hl,{'BC full support','Pressure surface 1','Pressure surface 2','BC force'});

axisGeom(gca,fontSize);
camlight headlight; 
gdrawnow; 

%% Split element regions (central plate and sides);

XE=patchCentre(E,V(:,1));

logicElementsInner=abs(XE)<eps;

E1=E(~logicElementsInner,:); %Other elements
E2=E(logicElementsInner,:); %Inner element layer
E=[E1;E2];

[F1]=element2patch(E1);
[F2]=element2patch(E2);

%%

cFigure; hold on; 
gpatch(F1,V,'bw','k',0.5);
gpatch(F2,V,'rw','k',0.5);
axisGeom; 
colormap(turbo(250)); icolorbar; 
camlight headlight; 
gdrawnow; 

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
febio_spec.Step.step{3}.Control=stepStruct.Control;
febio_spec.Step.step{3}.ATTR.id=3;
febio_spec.Step.step{4}.Control=stepStruct.Control;
febio_spec.Step.step{4}.ATTR.id=4;

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

materialName2='Material2';
febio_spec.Material.material{2}.ATTR.name=materialName2;
febio_spec.Material.material{2}.ATTR.type='Ogden';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.c1=c2;
febio_spec.Material.material{2}.m1=m2;
febio_spec.Material.material{2}.c2=c2;
febio_spec.Material.material{2}.m2=-m2;
febio_spec.Material.material{2}.k=k2;

materialName3='Material3';
febio_spec.Material.material{3}.ATTR.name=materialName3;
febio_spec.Material.material{3}.ATTR.type='rigid body';
febio_spec.Material.material{3}.ATTR.id=3;
febio_spec.Material.material{3}.density=1;
febio_spec.Material.material{3}.center_of_mass=center_of_mass;

%Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type 
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E1,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E1; %The element matrix

partName2='Part2';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='hex8'; %Element type 
febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=E2; %The element matrix

partName3='Part3';
febio_spec.Mesh.Elements{3}.ATTR.name=partName3; %Name of this part
febio_spec.Mesh.Elements{3}.ATTR.type='quad4'; %Element type 
febio_spec.Mesh.Elements{3}.elem.ATTR.id=size(E1,1)+size(E2,1)+(1:1:size(F_top,1))'; %Element id's
febio_spec.Mesh.Elements{3}.elem.VAL=F_top; %The element matrix

% -> Surfaces
surfaceName1='LoadedSurface1';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(F_pressure1,1))';
febio_spec.Mesh.Surface{1}.quad4.VAL=F_pressure1;

surfaceName2='LoadedSurface2';
febio_spec.Mesh.Surface{2}.ATTR.name=surfaceName2;
febio_spec.Mesh.Surface{2}.quad4.ATTR.id=(1:1:size(F_pressure2,1))';
febio_spec.Mesh.Surface{2}.quad4.VAL=F_pressure2;

% -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcSupportList);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat=materialName1;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.name=partName2;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat=materialName2;
febio_spec.MeshDomains.ShellDomain.ATTR.name=partName3;
febio_spec.MeshDomains.ShellDomain.ATTR.mat=materialName3;

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_xyz';
febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.x_dof=1;
febio_spec.Boundary.bc{1}.y_dof=1;
febio_spec.Boundary.bc{1}.z_dof=1;

% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_xyz';
febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.x_dof=1;
febio_spec.Boundary.bc{1}.y_dof=1;
febio_spec.Boundary.bc{1}.z_dof=1;
%Loads section
% -> Surface load
febio_spec.Loads.surface_load{1}.ATTR.type='pressure';
febio_spec.Loads.surface_load{1}.ATTR.surface=surfaceName1;
febio_spec.Loads.surface_load{1}.pressure.ATTR.lc=1;
febio_spec.Loads.surface_load{1}.pressure.VAL=appliedPressure1;
febio_spec.Loads.surface_load{1}.symmetric_stiffness=1;

febio_spec.Loads.surface_load{2}.ATTR.type='pressure';
febio_spec.Loads.surface_load{2}.ATTR.surface=surfaceName2;
febio_spec.Loads.surface_load{2}.pressure.ATTR.lc=1;
febio_spec.Loads.surface_load{2}.pressure.VAL=appliedPressure2;
febio_spec.Loads.surface_load{2}.symmetric_stiffness=1;

%Rigid section 
% ->Rigid body fix boundary conditions
febio_spec.Rigid.rigid_bc{1}.ATTR.name=['RigidFix_0',num2str(q)];
febio_spec.Rigid.rigid_bc{1}.ATTR.type='rigid_fixed';
febio_spec.Rigid.rigid_bc{1}.rb=3;
febio_spec.Rigid.rigid_bc{1}.Rx_dof=0;
febio_spec.Rigid.rigid_bc{1}.Ry_dof=1;
febio_spec.Rigid.rigid_bc{1}.Rz_dof=0;
febio_spec.Rigid.rigid_bc{1}.Ru_dof=0;
febio_spec.Rigid.rigid_bc{1}.Rv_dof=0;
febio_spec.Rigid.rigid_bc{1}.Rw_dof=0;

% ->Rigid body prescribe boundary conditions
s=[1 -1 -1];
for q=2:1:4
    febio_spec.Step.step{q}.Rigid.rigid_bc{1}.ATTR.name=['RigidPrescribe_0',num2str(q)];
    febio_spec.Step.step{q}.Rigid.rigid_bc{1}.ATTR.type='rigid_displacement';
    febio_spec.Step.step{q}.Rigid.rigid_bc{1}.rb=3;
    febio_spec.Step.step{q}.Rigid.rigid_bc{1}.dof='x';
    febio_spec.Step.step{q}.Rigid.rigid_bc{1}.value.ATTR.lc=q;
    febio_spec.Step.step{q}.Rigid.rigid_bc{1}.value.VAL=s(q-1).*prescribedDisplacement_X;
    febio_spec.Step.step{q}.Rigid.rigid_bc{1}.relative=1;
end

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.name='LC_1';
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
%febio_spec.LoadData.load_controller{1}.extend='CONSTANT';
febio_spec.LoadData.load_controller{1}.points.pt.VAL=[0 0; 1 1; 2 1; 3 1; 4 1];

febio_spec.LoadData.load_controller{2}.ATTR.name='LC_2';
febio_spec.LoadData.load_controller{2}.ATTR.id=2;
febio_spec.LoadData.load_controller{2}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate='LINEAR';
%febio_spec.LoadData.load_controller{2}.extend='CONSTANT';
febio_spec.LoadData.load_controller{2}.points.pt.VAL=[0 0; 1 0; 2 1;];

febio_spec.LoadData.load_controller{3}.ATTR.name='LC_3';
febio_spec.LoadData.load_controller{3}.ATTR.id=3;
febio_spec.LoadData.load_controller{3}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{3}.interpolate='LINEAR';
%febio_spec.LoadData.load_controller{3}.extend='CONSTANT';
febio_spec.LoadData.load_controller{3}.points.pt.VAL=[0 0; 1 0; 2 0; 3 1];

febio_spec.LoadData.load_controller{4}.ATTR.name='LC_4';
febio_spec.LoadData.load_controller{4}.ATTR.id=4;
febio_spec.LoadData.load_controller{4}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{4}.interpolate='LINEAR';
%febio_spec.LoadData.load_controller{4}.extend='CONSTANT';
febio_spec.LoadData.load_controller{4}.points.pt.VAL=[0 0; 1 0; 2 0; 3 0; 4 1];

%Output section 
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='s1';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.element_data{1}.VAL=1:size(E,1);

febio_spec.Output.logfile.rigid_body_data{1}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.rigid_body_data{1}.ATTR.data='Fx';
febio_spec.Output.logfile.rigid_body_data{1}.ATTR.delim=',';

%% Quick viewing of the FEBio input file structure
% The |febView| function can be used to view the xml structure in a MATLAB
% figure window. 

%%
% |febView(febio_spec); %Viewing the febio file|

%% Exporting the FEBio input file
% Exporting the febio_spec structure to an FEBio input file is done using
% the |febioStruct2xml| function. 

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
% febView(febioFebFileName); 

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
    % Importing nodal displacements from a log file
    dataStructForce=importFEBio_logfile(fullfile(savePath,febioLogFileName_force),0,1);
    timeData=dataStructForce.time(:);
    forceData=dataStructForce.data(:);
    
    logicPush=timeData>=1 & timeData<=2;
    timePush=timeData(logicPush);
    timePush=timePush-1;
    forcePush=forceData(logicPush);
    dispPush=timePush.*prescribedDisplacement_X;
    
    logicPull=timeData>=3 & timeData<=4;
    timePull=timeData(logicPull);
    timePull=timePull-3;
    forcePull=forceData(logicPull);
    dispPull=timePull.*-prescribedDisplacement_X;
    
    Ux=[flipud(dispPull(2:end)); dispPush(2:end)];
    Fx=[flipud(forcePull(2:end)); forcePush(2:end)];
    
    dF=diff(Fx);
    du=diff(Ux);
    S_diff=dF./du;
    u_diff=Ux(1:end-1)+du/2;
            
    ui=linspace(-prescribedDisplacement_X,prescribedDisplacement_X,100);
    S=interp1(u_diff,S_diff,ui,'linear','extrap');

    %%
    
    cFigure; hold on; 
    xlabel('U_x [mm]'); ylabel('F_x [N]');
    hp1=plot(dispPull,forcePull,'b.-','LineWidth',3,'MarkerSize',25);
    hp2=plot(dispPush,forcePush,'r.-','LineWidth',3,'MarkerSize',25);
    legend([hp1 hp2],{'Pull','Push'},'Location','NorthEastOutside')
    set(gca,'FontSize',fontSize);
    axis tight; axis square; grid on; box on;
    drawnow
    
    %%
    
    cFigure; hold on; 
    xlabel('U_x [mm]'); ylabel('S_x [N/mm]');
    hp1=plot(u_diff,S_diff,'k.','MarkerSize',50);
    hp2=plot(ui,S,'k-','LineWidth',3);
    legend([hp1 hp2],{'FEA','Interpolated'},'Location','NorthEastOutside')
    set(gca,'FontSize',fontSize);
    axis tight; axis square; grid on; box on;
    drawnow
    
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    DN_magnitude=sqrt(sum(N_disp_mat(:,:,end).^2,2)); %Current displacement magnitude
        
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('Displacement magnitude [mm]','Interpreter','Latex')
    hp=gpatch(Fb,V_DEF(:,:,end),DN_magnitude,'k',1); %Add graphics object to animate
%     hp.Marker='.';
%     hp.MarkerSize=markerSize2;
    hp.FaceColor='interp';
    
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
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),DN_magnitude}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;
            
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
%     hp.Marker='.';
%     hp.MarkerSize=markerSize2;
    hp.FaceColor='interp';
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([min(E_stress_mat(:)) max(E_stress_mat(:))]);    
    axis(axisLim(V_DEF)); %Set axis limits statically    
    camlight headlight;        
        
    % Set up animation features
    animStruct.Time=timeVec; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        
        [CV]=faceToVertexMeasure(E,V,E_stress_mat(:,:,qt));
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
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
