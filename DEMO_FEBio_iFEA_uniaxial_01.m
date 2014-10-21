%% DEMO_FEBio_iFEA_uniaxial_01
% Below is a demonstration for: 
% 1) Inverse FEA based material parameter optimisation 

%%

clear; close all; clc; 

%%
% Plot settings
figColor='w'; figColorDef='white';
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

modelName=fullfile(savePath,'iFEA_tempModel');

%Specifying dimensions and number of elements
sampleWidth=10;
sampleThickness=10; 
sampleHeight=10;
pointSpacings=2*ones(1,3);
initialArea=sampleWidth*sampleThickness;

numElementsWidth=round(sampleWidth/pointSpacings(1));
numElementsThickness=round(sampleThickness/pointSpacings(2));
numElementsHeight=round(sampleHeight/pointSpacings(3));

stretchLoad=0.7;
displacementMagnitude=[0 0 (stretchLoad*sampleHeight)-sampleHeight];

%True material parameter set
k_factor=1e4;
c1_true=1e-3; 
m1_true=12;
k_true=c1_true*k_factor; 

%Initial material parameter set
c1_ini=c1_true./2; 
m1_ini=m1_true+7;
k_ini=c1_ini*k_factor; 
P=[c1_ini m1_ini];

%% SIMULATE EXPERIMENTAL DATA

%Basic set
stress_cauchy_exp=[0;-0.0422226256000000;-0.0811346871800000;-0.119872916800000;-0.161466624000000;-0.209098742000000;-0.266409832800000;-0.337879334400000;-0.429344276800000;-0.548728823600000;-0.707119980000000];
stretch_exp=[1;0.970000000000000;0.940000000000000;0.910000000000000;0.880000000000000;0.850000000000000;0.820000000000000;0.790000000000000;0.760000000000000;0.730000000000000;0.700000000000000];

%Interpolate to higher sampling
n=100;
stretch_exp_n=linspace(1,stretchLoad,n);
stress_cauchy_exp_n = interp1(stretch_exp,stress_cauchy_exp,stretch_exp_n,'pchip');

%Override variables
stress_cauchy_exp=stress_cauchy_exp_n;
stretch_exp=stretch_exp_n;

%Add noise
stdNoise=0.01; %Standard deviation in units of stress
stress_cauchy_exp_n=stress_cauchy_exp_n+stdNoise.*randn(size(stress_cauchy_exp_n));

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
hf=figuremax(figColor,figColorDef);
title('Model surfaces','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',faceBoundaryMarker,'FaceAlpha',faceAlpha2,'lineWidth',edgeWidth,'edgeColor',edgeColor);

colormap(jet(6)); colorbar; 
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
drawnow; 

%% DEFINE BC's

%Define supported node sets
logicFace=faceBoundaryMarker==1;
Fr=Fb(logicFace,:);
bcSupportList_X=unique(Fr(:));

logicFace=faceBoundaryMarker==3;
Fr=Fb(logicFace,:);
bcSupportList_Y=unique(Fr(:));

logicFace=faceBoundaryMarker==5;
Fr=Fb(logicFace,:);
bcSupportList_Z=unique(Fr(:));

%Prescribed displacement nodes
logicPrescribe=faceBoundaryMarker==6;
Fr=Fb(logicPrescribe,:);
bcPrescribeList=unique(Fr(:));
bcPrescribeMagnitudes=displacementMagnitude(ones(1,numel(bcPrescribeList)),:);

%%
% Visualize BC's
hf=figuremax(figColor,figColorDef);
title('Complete model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',Fb,'Vertices',V,'FaceColor','flat','CData',faceBoundaryMarker,'FaceAlpha',faceAlpha2,'lineWidth',edgeWidth,'edgeColor',edgeColor);
plotV(V(bcSupportList_X,:),'r.','MarkerSize',markerSize);
plotV(V(bcSupportList_Y,:),'g.','MarkerSize',markerSize);
plotV(V(bcSupportList_Z,:),'b.','MarkerSize',markerSize);
plotV(V(bcPrescribeList,:),'k.','MarkerSize',markerSize);
set(gca,'FontSize',fontSize);

colormap(jet(6)); colorbar; 
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;
drawnow; 

%% CONSTRUCTING FEB MODEL

FEB_struct.febio_spec.version='2.0';

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

%Material 1 uncoupled hyperelastic
FEB_struct.Materials{1}.Type='Ogden';
FEB_struct.Materials{1}.Name='Block_mat';
FEB_struct.Materials{1}.Properties={'c1','m1','k'};
FEB_struct.Materials{1}.Values={c1_ini,m1_ini,k_ini};

%Control section
FEB_struct.Control.AnalysisType='static';
FEB_struct.Control.Properties={'time_steps','step_size',...
    'max_refs','max_ups',...
    'dtol','etol','rtol','lstol'};
FEB_struct.Control.Values={20,0.05,25,0,0.001,0.01,0,0.9};
FEB_struct.Control.TimeStepperProperties={'dtmin','dtmax','max_retries','opt_iter','aggressiveness'};
FEB_struct.Control.TimeStepperValues={1e-4,0.05,5,10,1};

%Defining node sets
FEB_struct.Geometry.NodeSet{1}.Set=bcSupportList_X;
FEB_struct.Geometry.NodeSet{1}.Name='bcSupportList_X';
FEB_struct.Geometry.NodeSet{2}.Set=bcSupportList_Y;
FEB_struct.Geometry.NodeSet{2}.Name='bcSupportList_Y';
FEB_struct.Geometry.NodeSet{3}.Set=bcSupportList_Z;
FEB_struct.Geometry.NodeSet{3}.Name='bcSupportList_Z';
% FEB_struct.Geometry.NodeSet{4}.Set=bcPrescribeList;
% FEB_struct.Geometry.NodeSet{4}.Name='bcPrescribeList';

%Adding BC information
FEB_struct.Boundary.Fix{1}.bc='x';
FEB_struct.Boundary.Fix{1}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
FEB_struct.Boundary.Fix{2}.bc='y';
FEB_struct.Boundary.Fix{2}.SetName=FEB_struct.Geometry.NodeSet{2}.Name;
FEB_struct.Boundary.Fix{3}.bc='z';
FEB_struct.Boundary.Fix{3}.SetName=FEB_struct.Geometry.NodeSet{3}.Name;

%Prescribed BC's
FEB_struct.Boundary.Prescribe{1}.Set=bcPrescribeList;
FEB_struct.Boundary.Prescribe{1}.bc='z';
FEB_struct.Boundary.Prescribe{1}.lc=1;
FEB_struct.Boundary.Prescribe{1}.nodeScale=displacementMagnitude(ones(numel(bcPrescribeList),1),3);
FEB_struct.Boundary.Prescribe{1}.Type='relative';

%Load curves
FEB_struct.LoadData.LoadCurves.id=1;
FEB_struct.LoadData.LoadCurves.type={'linear'};
FEB_struct.LoadData.LoadCurves.loadPoints={[0 0;1 1;]};

%Adding output requests
FEB_struct.Output.VarTypes={'displacement','stress','relative volume','shell thickness'};

%Specify log file output
run_disp_output_name=[FEB_struct.run_filename(1:end-4),'_node_out.txt'];
run_force_output_name=[FEB_struct.run_filename(1:end-4),'_force_out.txt'];
FEB_struct.run_output_names={run_disp_output_name,run_force_output_name};
FEB_struct.output_types={'node_data','node_data'};
FEB_struct.data_types={'ux;uy;uz','Rx;Ry;Rz'};

%% SAVING .FEB FILE
FEB_struct.disp_opt=0; %Display waitbars option
febStruct2febFile(FEB_struct);

%% RUNNING FEBIO JOB

FEBioRunStruct.FEBioPath='C:\Program Files\febio-2.1.0\bin\FEBio2.exe';
FEBioRunStruct.run_filename=FEB_struct.run_filename;
FEBioRunStruct.run_logname=FEB_struct.run_logname;
FEBioRunStruct.disp_on=1;
FEBioRunStruct.disp_log_on=1;
FEBioRunStruct.runMode='external';%'internal';
FEBioRunStruct.t_check=0.25; %Time for checking log file (dont set too small)
FEBioRunStruct.maxtpi=1e99; %Max analysis time
FEBioRunStruct.maxLogCheckTime=3; %Max log file checking time

[runFlag]=runMonitorFEBio(FEBioRunStruct);%START FEBio NOW!!!!!!!!

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

hf1=figuremax(figColor,figColorDef);
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
[time_mat, N_force_mat,~]=importFEBio_logfile(FEB_struct.run_output_names{2}); %Nodal displacements

FZ_set=N_force_mat(bcPrescribeList,end,:); %Final nodal displacements

%%

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

%Interpolate experiment onto simulated points
stress_cauchy_exp_sim = interp1(stretch_exp,stress_cauchy_exp,stretch_sim,'pchip');

%%

hf1=figuremax(figColor,figColorDef);
title('Stretch stress curves, initial differences between simulation and experiment','FontSize',fontSize);
xlabel('\lambda Stretch [.]','FontSize',fontSize); ylabel('\sigma Cauchy stress [kPa]','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;

H(1)=plot(stretch_sim,stress_cauchy_sim,'r.-','lineWidth',lineWidth,'markerSize',markerSize);
H(2)=plot(stretch_exp,stress_cauchy_exp,'g-','lineWidth',lineWidth);
H(3)=plot(stretch_sim,stress_cauchy_exp_sim,'k.','markerSize',markerSize);

legend(H,{'Simulation','Experiment','Experiment@simulation'},'Location','northwest');
view(2); axis tight;  grid on;
set(gca,'FontSize',fontSize);
drawnow;

%% CREATE STRUCTURES FOR OPTIMISATION 

mat_struct.id=1;
mat_struct.par_names={'c1','m1','k'};
mat_struct.par_values={c1_ini m1_ini k_ini}; 

% docNode=set_mat_par_FEBIO(FEB_struct.run_filename,FEB_struct.run_filename,{mat_struct});

FEBioRunStruct.disp_on=0; 
FEBioRunStruct.disp_log_on=0; 

%What should be known to the objective function:
objectiveStruct.bcPrescribeList=bcPrescribeList;
objectiveStruct.stretch_exp=stretch_exp;
objectiveStruct.stress_cauchy_exp=stress_cauchy_exp;
objectiveStruct.FEBioRunStruct=FEBioRunStruct;
objectiveStruct.FEB_struct=FEB_struct;
objectiveStruct.mat_struct=mat_struct;
objectiveStruct.k_factor=k_factor;
objectiveStruct.P_ini=P;
objectiveStruct.initialArea=initialArea;
objectiveStruct.sampleHeight=sampleHeight;

%Optimisation settings
maxNumberIterations=50; %Maximum number of iterations
maxNumberFunctionEvaluations=maxNumberIterations; %Maximum number of function evaluations, N.B. multiple evaluations are used per iteration
functionTolerance=1e-3; %Tolerance on objective function value
parameterTolerance=1e-3; %Tolerance on parameter variation
displayTypeIterations='iter';
OPT_options = optimset('MaxFunEvals',maxNumberFunctionEvaluations,...
                       'MaxIter',maxNumberIterations,...
                       'TolFun',functionTolerance,...
                       'TolX',parameterTolerance,...
                       'Display',displayTypeIterations);
objectiveStruct.method=1; 

%File names of output files
objectiveStruct.run_output_names=FEB_struct.run_output_names;

%% PARAMETER BOUNDS

Ps.f=1;
Ps.t=0.9;
Ps.ub=[0.1  22];
Ps.lb=[0  2];
Ps.c=P;
objectiveStruct.Ps=Ps;

Pb=zeros(size(P));
for q=1:1:numel(P);
    Psn=objectiveStruct.Ps;
    Psn.ub=objectiveStruct.Ps.ub(q);
    Psn.lb=objectiveStruct.Ps.lb(q);
    Psn.c=objectiveStruct.Ps.c(q);
    Pb(q)=inv_parbound(P(q),Psn);
end

%% STARTING OPTIMISATION

switch objectiveStruct.method
    case 1 %fminsearch and Nelder-Mead
        [Pb_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= fminsearch(@(Pb) obj_DEMO_FEBio_iFEA_uniaxial_01(Pb,objectiveStruct),Pb,OPT_options);
    case 2 %lsqnonlin and Levenberg-Marquardt
        OPT_options = optimset(OPT_options,'Algorithm',{'levenberg-marquardt',.01}); %Specifically setting algorithm
        [Pb_opt,OPT_out.resnorm,OPT_out.residual]= lsqnonlin(@(Pb) obj_DEMO_FEBio_iFEA_uniaxial_01(Pb,objectiveStruct),Pb,[],[],OPT_options);    
end

[Fopt,OPT_stats_out]=obj_DEMO_FEBio_iFEA_uniaxial_01(Pb_opt,objectiveStruct);

type(fullfile(fileparts(filePath),'obj_DEMO_FEBio_iFEA_uniaxial_01'))

%%
P_opt=zeros(size(P));
for q=1:1:numel(P);
    Psn=objectiveStruct.Ps;
    Psn.ub=objectiveStruct.Ps.ub(q);
    Psn.lb=objectiveStruct.Ps.lb(q);
    Psn.c=objectiveStruct.Ps.c(q);
    P_opt(q)=parbound(Pb_opt(q),Psn);
end

disp_text=sprintf('%6.16e,',P_opt); disp_text=disp_text(1:end-1);
disp(['P_opt=',disp_text]);

%%

hf1=figuremax(figColor,figColorDef);
title('Stretch stress curves optimised','FontSize',fontSize);
xlabel('\lambda Stretch [.]','FontSize',fontSize); ylabel('\sigma Cauchy stress [kPa]','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;

H(1)=plot(OPT_stats_out.stretch_sim,OPT_stats_out.stress_cauchy_sim,'r.-','lineWidth',lineWidth,'markerSize',markerSize);
H(2)=plot(stretch_exp,stress_cauchy_exp,'g-','lineWidth',lineWidth);

legend(H,{'Simulation','Experiment'},'Location','northwest');
view(2); axis tight;  grid on;
set(gca,'FontSize',fontSize);
drawnow;

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
