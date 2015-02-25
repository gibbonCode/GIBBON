%% DEMO_FEBio_block_uniaxial_compression_dynamic_viscoelastic_ramp_hold
% Below is a demonstration for:
% 1) Building an FEBio model for uniaxial compression for a viscoelastic
% material
% 2) Running the model
% 3) Importing displacement and force results
% 4) Plotting results

%%

clear; close all; clc;

%%
% Plot settings
fontSize=20;
faceAlpha1=0.8;
faceAlpha2=1;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize=25;
lineWidth=3;

%%
% Control parameters

% path names
filePath=mfilename('fullpath');
savePath=fullfile(fileparts(filePath),'data','temp');

modelName=fullfile(savePath,'tempModel');

%Specifying dimensions and number of elements
sampleWidth=10;
sampleThickness=10;
sampleHeight=10;
pointSpacings=[1 1 1]*2;
initialArea=sampleWidth*sampleThickness;

numElementsWidth=round(sampleWidth/pointSpacings(1));
numElementsThickness=round(sampleThickness/pointSpacings(2));
numElementsHeight=round(sampleHeight/pointSpacings(3));

stretchLoad=0.7;
displacementMagnitude=[0 0 (stretchLoad*sampleHeight)-sampleHeight];

%Material parameter set
c1=1e-3; %ogden c1
m1=3; %ogden m1
k_factor=1000; %Bulk modulus factor
k=c1*k_factor; %The bulk modulus
g1=1; %Viscoelastic QLV proportional coefficient
t1=0.25; %Viscoelastic QLV time coefficient
d=1e-9; %Density (not required for static analysis)
t_total=5; %Total simulation time
t_load=0.5; %Time from start to max load
t_step_ini=1e-3; %Initial desired step size
t_step_max=0.05; %Maximum step size

multiStep=0; 

uncoupledLaw=1; %1=uncoupled, 2=coupled

%% CREATING MESHED BOX

%Create box 1
boxDim=[sampleWidth sampleThickness sampleHeight]; %Dimensions
boxEl=[numElementsWidth numElementsThickness numElementsHeight]; %Number of elements
[box1]=hexMeshBox(boxDim,boxEl);
E=box1.E;
V=box1.V;
Fb=box1.Fb;
faceBoundaryMarker=box1.faceBoundaryMarker;

X=V(:,1); Y=V(:,2); Z=V(:,3);
VE=[mean(X(E),2) mean(Y(E),2) mean(Z(E),2)];

elementMaterialIndices=ones(size(E,1),1);

%%

% Plotting boundary surfaces
hf=cFigure;
title('Model surfaces','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',faceBoundaryMarker,'FaceAlpha',faceAlpha2,'lineWidth',edgeWidth,'edgeColor',edgeColor);

colormap(jet(6)); colorbar;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
drawnow;

%% DEFINE BC's

%Define faces
logicFace=faceBoundaryMarker==1;
Fr=Fb(logicFace,:);
bcSupportList_X=unique(Fr(:));

logicFace=faceBoundaryMarker==3;
Fr=Fb(logicFace,:);
bcSupportList_Y=unique(Fr(:));

logicFace=faceBoundaryMarker==5;
Fr=Fb(logicFace,:);
bcSupportList_Z=unique(Fr(:));

%Define line support
bcSupportList_X_axis=bcSupportList_Y(ismember(bcSupportList_Y,bcSupportList_Z));
bcSupportList_Y_axis=bcSupportList_X(ismember(bcSupportList_X,bcSupportList_Z));

%Prescribed displacement nodes
logicPrescribe=faceBoundaryMarker==6;
Fr=Fb(logicPrescribe,:);
bcPrescribeList=unique(Fr(:));
bcPrescribeMagnitudes=displacementMagnitude(ones(1,numel(bcPrescribeList)),:);

%%
% Visualize BC's
hf=cFigure;
title('Model BCs','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',faceBoundaryMarker,'FaceAlpha',faceAlpha2,'lineWidth',edgeWidth,'edgeColor',edgeColor);

plotV(V(bcSupportList_Z,:),'b.','MarkerSize',markerSize);
plotV(V(bcPrescribeList,:),'k.','MarkerSize',markerSize);
plotV(V(bcSupportList_X_axis,:),'g.','MarkerSize',markerSize);
plotV(V(bcSupportList_Y_axis,:),'r.','MarkerSize',markerSize);

set(gca,'FontSize',fontSize);
colormap(jet(6)); colorbar;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
drawnow;

%% CONSTRUCTING FEB MODEL

FEB_struct.febio_spec.version='2.0';
FEB_struct.Module.Type='solid';

% Defining file names
FEB_struct.run_filename=[modelName,'.feb']; %FEB file name
FEB_struct.run_logname=[modelName,'.txt']; %FEBio log file name

%Geometry section
FEB_struct.Geometry.Nodes=V;
FEB_struct.Geometry.Elements={E}; %The element sets
FEB_struct.Geometry.ElementType={'hex8'}; %The element types
FEB_struct.Geometry.ElementMat={elementMaterialIndices};
FEB_struct.Geometry.ElementsPartName={'Block'};

%Material section
switch uncoupledLaw
    case 1
        FEB_struct.Materials{1}.Type='uncoupled viscoelastic';
        FEB_struct.Materials{1}.Name='Block_material';
        FEB_struct.Materials{1}.Properties={'g1','t1','elastic','density'};
        FEB_struct.Materials{1}.Values={g1,t1,[],d};
        FEB_struct.Materials{1}.PropAttrName=cell(1,numel(FEB_struct.Materials{1}.Properties));
        FEB_struct.Materials{1}.PropAttrName{3}='type';
        FEB_struct.Materials{1}.PropAttrVal{3}='Ogden';
        FEB_struct.Materials{1}.PropParName=cell(1,numel(FEB_struct.Materials{1}.Properties));
        FEB_struct.Materials{1}.PropParVal=cell(1,numel(FEB_struct.Materials{1}.Properties));
        FEB_struct.Materials{1}.PropParName{3}={'c1','m1','k','density'};
        FEB_struct.Materials{1}.PropParVal{3}={c1,m1,k,d};
    case 0
        FEB_struct.Materials{1}.Type='viscoelastic';
        FEB_struct.Materials{1}.Name='Block_material';
        FEB_struct.Materials{1}.Properties={'g1','t1','elastic','density'};
        FEB_struct.Materials{1}.Values={g1,t1,[],d};
        FEB_struct.Materials{1}.PropAttrName=cell(1,numel(FEB_struct.Materials{1}.Properties));
        FEB_struct.Materials{1}.PropAttrName{3}='type';
        FEB_struct.Materials{1}.PropAttrVal{3}='Ogden unconstrained';
        FEB_struct.Materials{1}.PropParName=cell(1,numel(FEB_struct.Materials{1}.Properties));
        FEB_struct.Materials{1}.PropParVal=cell(1,numel(FEB_struct.Materials{1}.Properties));
        FEB_struct.Materials{1}.PropParName{3}={'c1','m1','cp','density'};
        FEB_struct.Materials{1}.PropParVal{3}={c1,m1,k,d};
end

%Step specific control sections
n=round(t_total/t_step_ini);
t_step=t_total/n;
FEB_struct.Control.AnalysisType='static';
FEB_struct.Control.Properties={'time_steps','step_size',...
    'max_refs','max_ups',...
    'dtol','etol','rtol','lstol'};
FEB_struct.Control.Values={n,t_step,...
    15,0,...
    0.001,0.01,0,0.9};
FEB_struct.Control.TimeStepperProperties={'dtmin','dtmax','max_retries','opt_iter','aggressiveness'};
FEB_struct.Control.TimeStepperValues={t_step/100,t_step_max,5,10,0};

%Defining node sets
%Defining node sets
FEB_struct.Geometry.NodeSet{1}.Set=bcSupportList_Y_axis;
FEB_struct.Geometry.NodeSet{1}.Name='bcSupportList_Y_axis';
FEB_struct.Geometry.NodeSet{2}.Set=bcSupportList_X_axis;
FEB_struct.Geometry.NodeSet{2}.Name='bcSupportList_X_axis';
FEB_struct.Geometry.NodeSet{3}.Set=bcSupportList_Z;
FEB_struct.Geometry.NodeSet{3}.Name='bcSupportList_Z';
FEB_struct.Geometry.NodeSet{4}.Set=bcPrescribeList;
FEB_struct.Geometry.NodeSet{4}.Name='bcPrescribeList';

%Adding BC information
FEB_struct.Boundary.Fix{1}.bc='x';
FEB_struct.Boundary.Fix{1}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
FEB_struct.Boundary.Fix{2}.bc='y';
FEB_struct.Boundary.Fix{2}.SetName=FEB_struct.Geometry.NodeSet{2}.Name;
FEB_struct.Boundary.Fix{3}.bc='z';
FEB_struct.Boundary.Fix{3}.SetName=FEB_struct.Geometry.NodeSet{3}.Name;

FEB_struct.Boundary.Prescribe{1}.SetName=FEB_struct.Geometry.NodeSet{4}.Name;
FEB_struct.Boundary.Prescribe{1}.Scale=displacementMagnitude(3);
FEB_struct.Boundary.Prescribe{1}.bc='z';
FEB_struct.Boundary.Prescribe{1}.lc=1;
FEB_struct.Boundary.Prescribe{1}.Type='relative';

%Load curves
FEB_struct.LoadData.LoadCurves.id=1;
FEB_struct.LoadData.LoadCurves.type={'linear'};
FEB_struct.LoadData.LoadCurves.loadPoints={[0 0;t_load 1;t_total 1]};

%Adding output requests
FEB_struct.Output.VarTypes={'displacement','stress','relative volume','shell thickness'};

%Specify log file output
run_disp_output_name=[FEB_struct.run_filename(1:end-4),'_node_out.txt'];
run_force_output_name=[FEB_struct.run_filename(1:end-4),'_force_out.txt'];
FEB_struct.run_output_names={run_disp_output_name,run_force_output_name};
FEB_struct.output_types={'node_data','node_data'};
FEB_struct.data_types={'ux;uy;uz','Rx;Ry;Rz'};

%% SAVING .FEB FILE

FEB_struct.disp_opt=0; %Display waitbars
febStruct2febFile(FEB_struct);

%% RUNNING FEBIO JOB

FEBioRunStruct.FEBioPath='C:\Program Files\febio-2.2.2\bin\FEBio2.exe';
FEBioRunStruct.run_filename=FEB_struct.run_filename;
FEBioRunStruct.run_logname=FEB_struct.run_logname;
FEBioRunStruct.disp_on=1;
FEBioRunStruct.disp_log_on=1;
FEBioRunStruct.runMode='external';%'internal';
FEBioRunStruct.t_check=0.25; %Time for checking log file (dont set too small)
FEBioRunStruct.maxtpi=1e99; %Max analysis time
FEBioRunStruct.maxLogCheckTime=3; %Max log file checking time

[runFlag]=runMonitorFEBio(FEBioRunStruct);%START FEBio NOW!!!!!!!!

%%
if runFlag==1 %i.e. a succesful run
    
    %% IMPORTING NODAL DISPLACEMENT RESULTS
    % Importing nodal displacements from a log file
    [~, N_disp_mat,~]=importFEBio_logfile(FEB_struct.run_output_names{1}); %Nodal displacements
    
    DN=N_disp_mat(:,2:end,end); %Final nodal displacements
    
    %% CREATING NODE SET IN DEFORMED STATE
    V_def=V+DN;
    DN_magnitude=sqrt(sum(DN.^2,2));
    
    %%
    % Plotting the deformed model
    
    [CF]=vertexToFaceMeasure(Fb,DN_magnitude);
    
    hf1=cFigure;
    title('The deformed model','FontSize',fontSize);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
    
    hps=patch('Faces',Fb,'Vertices',V_def,'FaceColor','flat','CData',CF);
    
    view(3); axis tight;  axis equal;  grid on;
    colormap jet; colorbar;
    % camlight headlight;
    set(gca,'FontSize',fontSize);
    drawnow;
    
    %% IMPORTING NODAL FORCES
    % Importing nodal forces from a log file
    [time_mat, N_force_mat,~]=importFEBio_logfile(FEB_struct.run_output_names{2}); %Nodal forces
    time_mat=[0; time_mat(:)]; %Time
    
    %% DERIVING STRESS METRICS
    
    %Get Z forces
    FZ=sum(N_force_mat(bcPrescribeList,end,:),1);
    FZ=[0; FZ(:)]; %Mean top surface nodal forces
    
    %Derive applied stretch
    DZ_set=N_disp_mat(bcPrescribeList,end,:); %Final nodal displacements
    DZ_set=mean(DZ_set,1);
    stretch_sim=(DZ_set+sampleHeight)./sampleHeight;
    stretch_sim=[1; stretch_sim(:)];
    
    %Derive simulated Cauchy stress (alternatively import stress and take the mean)
    currentArea=initialArea./stretch_sim;
    stress_cauchy_sim=FZ./currentArea; %Cauchy stress
    stress_cauchy_sim=stress_cauchy_sim.*1e3; %Scale to kPa
    
    %%
    
    hf1=cFigure;
    title('Stress relaxation curve','FontSize',fontSize);
    xlabel('Time (s)','FontSize',fontSize); ylabel('\sigma Cauchy stress (kPa)','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
    
    plot(time_mat(:),stress_cauchy_sim(:),'r.-','lineWidth',lineWidth,'markerSize',markerSize);
    
    view(2); axis tight;  grid on;
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
