%% DEMO_febio_0089_iFEA_goh_skin_01
% Below is a demonstration for:
% 
% * Inverse FEA based material parameter optimisation for uniaxial
% compression

%% Keywords
%
% * febio_spec version 4.0
% * febio, FEBio
% * hexahedral elements, hex8
% * static, solid
% * hyperelastic, Ogden
% * displacement logfile
% * stress logfile

%%

clear; close all; clc;

%%
% Plot settings
fontSize=20;
faceAlpha1=0.8;
faceAlpha2=1;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize1=15;
markerSize2=40;
markerSize3=50;
lineWidth1=3; 
lineWidth2=3; 
cMap=viridis(20);
numPlotPointsGraphs = 100;

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress sigma_z
febioLogFileName_stretch=[febioFebFileNamePart,'_stretch_out.txt']; %Log file name for exporting stretch U_z

% Define data paths
loadPath_experimental = fullfile(defaultFolder,'data','Ni_Annaidh_2012');
dataName_1 = fullfile(loadPath_experimental,'Ni_Annaidh_stress_stretch_perp.mat');
dataName_2 = fullfile(loadPath_experimental,'Ni_Annaidh_stress_stretch_para.mat');

%Specifying dimensions and number of elements
sampleWidth=10;
sampleThickness=10; 
sampleHeight=10;
pointSpacings=10*ones(1,3);
initialArea=sampleWidth*sampleThickness;

numElementsWidth=round(sampleWidth/pointSpacings(1));
numElementsThickness=round(sampleThickness/pointSpacings(2));
numElementsHeight=round(sampleHeight/pointSpacings(3));

%Initial material parameter set
k_factor  = 1000;
c_ini     = 1.3021255296772; % Initial slope
k1_ini    = 124.3290999878474; % "Modulus of fibers"
k2_ini    = 0.0001368087512; %  "non-linearity" or J-shape
kappa_ini = 0.2985113701771; % [0,1/3], "distance between graphs"
k_ini     = c_ini*k_factor; % Bulk modulus
gamma     = 41;

evalMode = 1; % 1= test, 2=optimise

P = [c_ini k1_ini k2_ini kappa_ini]; % Parameter vector
% P = [c_ini*0.9 k1_ini*0.9 k2_ini*0.9 kappa_ini*0.9]; % Parameter vector
% P = [c_ini*1.1 k1_ini*1.1 k2_ini*1.1 kappa_ini*1.1]; % Parameter vector

LangerAngle_perp = 0/180*pi; % Angle of Langer line with loading direction
LangerAngle_para = 90/180*pi; % Angle of Langer line with loading direction 

LangerAngles = [LangerAngle_perp, LangerAngle_para];

% FEA control settings
numTimeSteps=20; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

runMode='external';% 'internal' or 'external'

%Optimisation settings
maxNumberIterations=100; %Maximum number of optimization iterations
maxNumberFunctionEvaluations=maxNumberIterations*10; %Maximum number of function evaluations, N.B. multiple evaluations are used per iteration
functionTolerance=1e-6; %Tolerance on objective function value
parameterTolerance=1e-6; %Tolerance on parameter variation
displayTypeIterations='iter';
optimisationMethod = 2; % 1= fminsearch, Nelder-Mead, 2=lsqnonlin, Levenberg-Marquart

%% LOAD EXPERIMENTAL DATA

data_perp = load(dataName_1);
data_para = load(dataName_2);

stretch_exp_perp_raw = data_perp.X1;
stress_exp_perp_raw = data_perp.Y1;

stretch_exp_para_raw = data_para.X1;
stress_exp_para_raw = data_para.Y1;

%%
% Set stretch levels for simulation
stretchLoad_perp = max(stretch_exp_perp_raw);
stretchLoad_para = max(stretch_exp_para_raw);

displacementMagnitude_perp=(stretchLoad_perp*sampleHeight)-sampleHeight;
displacementMagnitude_para=(stretchLoad_para*sampleHeight)-sampleHeight;

displacementMagnitudes = [displacementMagnitude_perp, displacementMagnitude_para];

%% 

stretch_plot_perp = linspace(1,stretchLoad_perp,numPlotPointsGraphs);
stretch_plot_para = linspace(1,stretchLoad_para,numPlotPointsGraphs);

stress_plot_perp = interp1(stretch_exp_perp_raw,stress_exp_perp_raw,stretch_plot_perp,'pchip','extrap'); %ppval(pp_perp,stretch_plot_perp);
stress_plot_para = interp1(stretch_exp_para_raw,stress_exp_para_raw,stretch_plot_para,'pchip','extrap'); %ppval(pp_para,stretch_plot_para);

%%

hf1 = cFigure; hold on;
title('Stretch stress curves optimised','FontSize',fontSize);
xlabel('\lambda Stretch [.]','FontSize',fontSize); ylabel('\sigma Cauchy stress [MPa]','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;

Hn(1) = plot(stretch_exp_perp_raw,stress_exp_perp_raw,'r.','MarkerSize',markerSize2);
Hn(2) = plot(stretch_plot_perp,stress_plot_perp,'r-','lineWidth',lineWidth1);

Hn(3) = plot(stretch_exp_para_raw,stress_exp_para_raw,'g.','MarkerSize',markerSize2);
Hn(4) = plot(stretch_plot_para,stress_plot_para,'g-','lineWidth',lineWidth1);

Hn(5) = plot(NaN,NaN,'ko','MarkerSize',markerSize1,'lineWidth',lineWidth2,'MarkerFaceColor','r');
Hn(6) = plot(NaN,NaN,'ko','MarkerSize',markerSize1,'lineWidth',lineWidth2,'MarkerFaceColor','g');

Hn(7) = plot(NaN,NaN,'rx','MarkerSize',markerSize1);
Hn(8) = plot(NaN,NaN,'gx','MarkerSize',markerSize1);

legend(Hn,{'Ni Annaidh - perp','interp - perp','Ni Annaidh - para','interp - para','FEA perp','FEA para'},'Location','northwest');
view(2); axis tight;  axis square; axis manual; box on; grid on; 
set(gca,'FontSize',fontSize);
drawnow;

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

cFigure; hold on;
title('Model surfaces','FontSize',fontSize);
gpatch(Fb,V,faceBoundaryMarker,'k',0.5);
colormap(gjet(6)); icolorbar; 
axisGeom(gca,fontSize);
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

%%
% Visualize BC's
cFigure; hold on;
title('Complete model','FontSize',fontSize);

gpatch(Fb,V,'kw','k',0.5);
plotV(V(bcSupportList_X,:),'r.','MarkerSize',markerSize1);
plotV(V(bcSupportList_Y,:),'g.','MarkerSize',markerSize1);
plotV(V(bcSupportList_Z,:),'b.','MarkerSize',markerSize1);
plotV(V(bcPrescribeList,:),'k.','MarkerSize',markerSize1);

axisGeom(gca,fontSize);
drawnow; 

%% DEFINE FIBRE DIRECTIONS

R = euler2DCM([0,-LangerAngle_perp,0]);
e1 = (R*[1 0 0]')';
e1_dir_perp = e1(ones(size(E,1),1),:);

e2 = [0 1 0];
e2_dir_perp = e2(ones(size(E,1),1),:);
e3_dir_perp = cross(e1_dir_perp,e2_dir_perp,2);

R = euler2DCM([0,-LangerAngle_para,0]);
e1 = (R*[1 0 0]')';
e1_dir_para = e1(ones(size(E,1),1),:);

e2 = [0 1 0];
e2_dir_para = e2(ones(size(E,1),1),:);
e3_dir_para = cross(e1_dir_para,e2_dir_para,2);

%% 
% Visualizing material directions

[VE]=patchCentre(E,V);

hf=cFigure;
subplot(1,2,1);  hold on;
title('Material directions - perp','FontSize',fontSize);

gpatch(Fb,V,'kw','none',0.25);
hf(1)=quiverVec(VE,e1_dir_perp,mean(pointSpacings),'r');
hf(2)=quiverVec(VE,e2_dir_perp,mean(pointSpacings)/2,'g');
hf(3)=quiverVec(VE,e3_dir_perp,mean(pointSpacings)/2,'b');

legend(hf,{'e1-direction','e2-direction','e3-direction'});
axisGeom(gca,fontSize);
camlight headlight; 

subplot(1,2,2);  hold on;
title('Material directions - para','FontSize',fontSize);

gpatch(Fb,V,'kw','none',0.25);
hf(1)=quiverVec(VE,e1_dir_para,mean(pointSpacings),'r');
hf(2)=quiverVec(VE,e2_dir_para,mean(pointSpacings)/2,'g');
hf(3)=quiverVec(VE,e3_dir_para,mean(pointSpacings)/2,'b');

legend(hf,{'e1-direction (fiber)','e2-direction','e3-direction'});
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

%Control section
febio_spec.Control.analysis='STATIC';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.solver.qn_method.max_ups=max_ups;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;

%Material section
materialName1='Material1';
febio_spec.Material.material{1}.ATTR.name=materialName1;
febio_spec.Material.material{1}.ATTR.type='HGO unconstrained';
% febio_spec.Material.material{1}.ATTR.type='Holzapfel-Gasser-Ogden';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c=c_ini;
febio_spec.Material.material{1}.k1=k1_ini;
febio_spec.Material.material{1}.k2=k2_ini;
febio_spec.Material.material{1}.gamma=gamma;
febio_spec.Material.material{1}.kappa=kappa_ini;
febio_spec.Material.material{1}.k=k_ini;

% Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix
 
% -> NodeSets
nodeSetName1='bcSupportList_X';
nodeSetName2='bcSupportList_Y';
nodeSetName3='bcSupportList_Z';
nodeSetName4='bcPrescribeList';

febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcSupportList_X);

febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.VAL=mrow(bcSupportList_Y);

febio_spec.Mesh.NodeSet{3}.ATTR.name=nodeSetName3;
febio_spec.Mesh.NodeSet{3}.VAL=mrow(bcSupportList_Z);

febio_spec.Mesh.NodeSet{4}.ATTR.name=nodeSetName4;
febio_spec.Mesh.NodeSet{4}.VAL=mrow(bcPrescribeList);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

%MeshData section
% -> ElementData
febio_spec.MeshData.ElementData{1}.ATTR.elem_set=partName1;
febio_spec.MeshData.ElementData{1}.ATTR.type='mat_axis';

for q=1:1:size(E,1)
    febio_spec.MeshData.ElementData{1}.elem{q}.ATTR.lid=q;
    febio_spec.MeshData.ElementData{1}.elem{q}.a=e1_dir_perp(q,:);
    febio_spec.MeshData.ElementData{1}.elem{q}.d=e2_dir_perp(q,:);
end

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_x';
febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.x_dof=1;
febio_spec.Boundary.bc{1}.y_dof=0;
febio_spec.Boundary.bc{1}.z_dof=0;

febio_spec.Boundary.bc{2}.ATTR.name='zero_displacement_y';
febio_spec.Boundary.bc{2}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{2}.x_dof=0;
febio_spec.Boundary.bc{2}.y_dof=1;
febio_spec.Boundary.bc{2}.z_dof=0;

febio_spec.Boundary.bc{3}.ATTR.name='zero_displacement_z';
febio_spec.Boundary.bc{3}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName3;
febio_spec.Boundary.bc{3}.x_dof=0;
febio_spec.Boundary.bc{3}.y_dof=0;
febio_spec.Boundary.bc{3}.z_dof=1;

febio_spec.Boundary.bc{4}.ATTR.name='prescibed_displacement_z';
febio_spec.Boundary.bc{4}.ATTR.type='prescribed displacement';
febio_spec.Boundary.bc{4}.ATTR.node_set=nodeSetName4;
febio_spec.Boundary.bc{4}.dof='z';
febio_spec.Boundary.bc{4}.value.ATTR.lc=1;
febio_spec.Boundary.bc{4}.value.VAL=displacementMagnitude_perp;
febio_spec.Boundary.bc{4}.relative=0;

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

febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='sz';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{2}.ATTR.file=febioLogFileName_stretch;
febio_spec.Output.logfile.element_data{2}.ATTR.data='Uz';
febio_spec.Output.logfile.element_data{2}.ATTR.delim=',';

% Plotfile section
febio_spec.Output.plotfile.compression=0;

%% Creating febio analysis structure 

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=0; %Display information on the command window
febioAnalysis.runMode=runMode;
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

%% Create structures for optimization 

%What should be known to the objective function:
objectiveStruct.stretch_exp_perp_raw=stretch_exp_perp_raw;
objectiveStruct.stress_exp_perp_raw=stress_exp_perp_raw;
objectiveStruct.stretch_exp_para_raw=stretch_exp_para_raw;
objectiveStruct.stress_exp_para_raw=stress_exp_para_raw;

objectiveStruct.febioAnalysis=febioAnalysis;
objectiveStruct.febio_spec=febio_spec;
objectiveStruct.febioFebFileName=febioFebFileName;
objectiveStruct.febioLogFileName_disp = fullfile(savePath,febioLogFileName_disp);
objectiveStruct.febioLogFileName_stress = fullfile(savePath,febioLogFileName_stress);
objectiveStruct.febioLogFileName_stretch = fullfile(savePath,febioLogFileName_stretch);

objectiveStruct.parNormFactors=P; %This will normalize the parameters to ones(size(P))
objectiveStruct.Pb_struct.xx_c=P; %Parameter constraining centre
objectiveStruct.Pb_struct.xxlim=[P(1)./1000 P(1).*1000;... % c
                                 P(2)./1000 P(2).*1000;... % k1
                                 0 P(3).*1000;... % k2
                                 0 1/3;... % kappa
                                    ]; %Parameter bounds
objectiveStruct.k_factor=k_factor;

objectiveStruct.LangerAngles = LangerAngles;
objectiveStruct.displacementMagnitudes = displacementMagnitudes;

objectiveStruct.hf = hf1;
objectiveStruct.h=[Hn(5),Hn(6),Hn(7),Hn(8)];

objectiveStruct.method = optimisationMethod; 

Pn=P./objectiveStruct.parNormFactors;

switch evalMode
    case 1
        [errorVal,simData] = objectiveFunctionIFEA(Pn, objectiveStruct);
    case 2
    %% start optimization
    switch optimisationMethod
        case 1 %fminsearch and Nelder-Mead
            OPT_options=optimset('fminsearch'); % 'Nelder-Mead simplex direct search'
            OPT_options = optimset(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
                'MaxIter',maxNumberIterations,...
                'TolFun',functionTolerance,...
                'TolX',parameterTolerance,...
                'Display',displayTypeIterations,...
                'FinDiffRelStep',1e-2,...
                'DiffMaxChange',0.5);
            [Pn_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= fminsearch(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),Pn,OPT_options);
        case 2 %lsqnonlin and Levenberg-Marquardt
            OPT_options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
            OPT_options = optimoptions(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
                'MaxIter',maxNumberIterations,...
                'TolFun',functionTolerance,...
                'TolX',parameterTolerance,...
                'Display',displayTypeIterations,...
                'FinDiffRelStep',1e-2,...
                'DiffMaxChange',0.5);
            [Pn_opt,OPT_out.resnorm,OPT_out.residual]= lsqnonlin(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),Pn,[],[],OPT_options);    
    end
    
    %% Unnormalize and constrain parameters
    
    [errorVal,simData] = objectiveFunctionIFEA(Pn_opt, objectiveStruct);
    
    P_opt=Pn_opt.*objectiveStruct.parNormFactors; %Scale back, undo normalization
    
    %Constraining parameters
    for q=1:1:numel(P_opt)
        [P_opt(q)]=boxconstrain(P_opt(q),objectiveStruct.Pb_struct.xxlim(q,1),objectiveStruct.Pb_struct.xxlim(q,2),objectiveStruct.Pb_struct.xx_c(q));    
    end
    
    disp_text=sprintf('%6.16e,',P_opt); disp_text=disp_text(1:end-1);
    disp(['P_opt=',disp_text]);

end

%% Import FEBio results 

% Displacements
N_disp_mat1=simData(1).dataStruct_disp.data; %Displacement
timeVec1=simData(1).dataStruct_disp.time; %Time
V_DEF1=N_disp_mat1+repmat(V,[1 1 size(N_disp_mat1,3)]); %Deformed coordinate set
DN_magnitude1=sqrt(sum(N_disp_mat1(:,:,end).^2,2)); %Current displacement magnitude

N_disp_mat2=simData(2).dataStruct_disp.data; %Displacement
timeVec2=simData(2).dataStruct_disp.time; %Time
V_DEF2=N_disp_mat2+repmat(V,[1 1 size(N_disp_mat2,3)]); %Deformed coordinate set
DN_magnitude2=sqrt(sum(N_disp_mat2(:,:,end).^2,2)); %Current displacement magnitude

%% 
% Plotting the simulated results using |anim8| to visualize and animate
% deformations 
   
% Create basic view and store graphics handle to initiate animation
hf=cFigure; %Open figure  
subplot(1,2,1);
title('Displ. magnitude [mm]','Interpreter','Latex')
hp1=gpatch(Fb,V_DEF1(:,:,end),DN_magnitude1,'k',1,2); %Add graphics object to animate
hp1.Marker='.'; hp1.MarkerSize=markerSize2; hp1.FaceColor='interp';
gpatch(Fb,V,0.5*ones(1,3),'none',0.25); %A static graphics object

axisGeom(gca,fontSize); 
colormap(cMap); colorbar;
caxis([0 max([max(DN_magnitude1) max(DN_magnitude2)])]); caxis manual;   
axis(axisLim([V_DEF1; V_DEF2])); %Set axis limits statically    
view(140,30);
camlight headlight;        
    
subplot(1,2,2);
title('Displ. magnitude [mm]','Interpreter','Latex')
hp2=gpatch(Fb,V_DEF2(:,:,end),DN_magnitude2,'k',1,2); %Add graphics object to animate
hp2.Marker='.'; hp2.MarkerSize=markerSize2; hp2.FaceColor='interp';
gpatch(Fb,V,0.5*ones(1,3),'none',0.25); %A static graphics object

axisGeom(gca,fontSize); 
colormap(cMap); colorbar;
caxis([0 max([max(DN_magnitude1) max(DN_magnitude2)])]); caxis manual;   
axis(axisLim([V_DEF1; V_DEF2])); %Set axis limits statically    
view(140,30);
camlight headlight;        

% Set up animation features
animStruct.Time=timeVec1; %The time vector    
for qt=1:1:size(N_disp_mat1,3) %Loop over time increments        
    
    DN_magnitude1=sqrt(sum(N_disp_mat1(:,:,qt).^2,2)); %Current displacement magnitude
    DN_magnitude2=sqrt(sum(N_disp_mat2(:,:,qt).^2,2)); %Current displacement magnitude        

    %Set entries in animation structure
    animStruct.Handles{qt}=[hp1 hp1 hp2 hp2]; %Handles of objects to animate
    animStruct.Props{qt}={'Vertices','CData','Vertices','CData'}; %Properties of objects to animate
    animStruct.Set{qt}={V_DEF1(:,:,qt),DN_magnitude1,V_DEF2(:,:,qt),DN_magnitude2}; %Property values for to set in order to animate
end        
anim8(hf,animStruct); %Initiate animation feature    
drawnow;
        
%%

%Access data
E_stress_mat1=simData(1).dataStruct_stress.data;
E_stretch_mat1=simData(1).dataStruct_stretch.data;
[CV1]=faceToVertexMeasure(E,V,E_stress_mat1(:,:,end));

E_stress_mat2=simData(2).dataStruct_stress.data;
E_stretch_mat2=simData(2).dataStruct_stretch.data;
[CV2]=faceToVertexMeasure(E,V,E_stress_mat2(:,:,end));

%% 
% Plotting the simulated results using |anim8| to visualize and animate
% deformations 

% Create basic view and store graphics handle to initiate animation
hf=cFigure; %Open figure  /usr/local/MATLAB/R2020a/bin/glnxa64/jcef_helper: symbol lookup error: /lib/x86_64-linux-gnu/libpango-1.0.so.0: undefined symbol: g_ptr_array_copy
subplot(1,2,1);
title('$\sigma_{zz}$ [MPa]','Interpreter','Latex')
hp1=gpatch(Fb,V_DEF1(:,:,end),CV1,'k',1,2); %Add graphics object to animate
hp1.Marker='.'; hp1.MarkerSize=markerSize2; hp1.FaceColor='interp';
gpatch(Fb,V,0.5*ones(1,3),'none',0.25); %A static graphics object

axisGeom(gca,fontSize); 
colormap(gca,cMap); colorbar;
caxis([min(E_stress_mat1(:)) max(E_stress_mat1(:))]);    
axis(axisLim([V_DEF1; V_DEF2])); %Set axis limits statically    
view(140,30);
camlight headlight;        
    
subplot(1,2,2);
title('$\sigma_{zz}$ [MPa]','Interpreter','Latex')
hp2=gpatch(Fb,V_DEF2(:,:,end),CV2,'k',1,2); %Add graphics object to animate
hp2.Marker='.'; hp2.MarkerSize=markerSize2; hp2.FaceColor='interp';
gpatch(Fb,V,0.5*ones(1,3),'none',0.25); %A static graphics object

axisGeom(gca,fontSize); 
colormap(gca,cMap); colorbar;
caxis([min(E_stress_mat2(:)) max(E_stress_mat2(:))]);    
axis(axisLim([V_DEF1; V_DEF2])); %Set axis limits statically    
view(140,30);
camlight headlight;        

% Set up animation features
animStruct.Time=timeVec1; %The time vector    
for qt=1:1:size(N_disp_mat1,3) %Loop over time increments        
    
    [CV1]=faceToVertexMeasure(E,V,E_stress_mat1(:,:,qt));
    [CV2]=faceToVertexMeasure(E,V,E_stress_mat2(:,:,qt));
    
    %Set entries in animation structure
    animStruct.Handles{qt}=[hp1 hp1 hp2 hp2]; %Handles of objects to animate
    animStruct.Props{qt}={'Vertices','CData','Vertices','CData'}; %Properties of objects to animate
    animStruct.Set{qt}={V_DEF1(:,:,qt),CV1,V_DEF2(:,:,qt),CV2}; %Property values for to set in order to animate
end        
anim8(hf,animStruct); %Initiate animation feature    
drawnow;

%%

function [errorVal,simData] = objectiveFunctionIFEA(Pn, objectiveStruct)

%% Access input structure
stretch_exp_perp_raw = objectiveStruct.stretch_exp_perp_raw;
stress_exp_perp_raw  = objectiveStruct.stress_exp_perp_raw;
stretch_exp_para_raw = objectiveStruct.stretch_exp_para_raw;
stress_exp_para_raw = objectiveStruct.stress_exp_para_raw;

febioAnalysis = objectiveStruct.febioAnalysis;
febio_spec = objectiveStruct.febio_spec;
febioFebFileName = objectiveStruct.febioFebFileName;
febioLogFileName_disp = objectiveStruct.febioLogFileName_disp;
febioLogFileName_stress = objectiveStruct.febioLogFileName_stress;
febioLogFileName_stretch = objectiveStruct.febioLogFileName_stretch;

parNormFactors = objectiveStruct.parNormFactors; 
Pb_struct = objectiveStruct.Pb_struct;
k_factor = objectiveStruct.k_factor;

LangerAngles = objectiveStruct.LangerAngles;
displacementMagnitudes = objectiveStruct.displacementMagnitudes;

hf = objectiveStruct.hf;
h = objectiveStruct.h;

%% Set/get material parameters
% Unnormalize and constrain parameters
P = Pn.*parNormFactors; %Scale back, undo normalization

%Constraining parameters
for q=1:1:numel(P)
    [P(q)]=boxconstrain(P(q),Pb_struct.xxlim(q,1),Pb_struct.xxlim(q,2),objectiveStruct.Pb_struct.xx_c(q));    
end

disp(['Trying   : ',sprintf(repmat('%6.8e ',[1,numel(P)]),P)]);

c = P(1); 
k1 = P(2);
k2 = P(3);
kappa = P(4);
k = c*k_factor;

% Set material parameters 
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c=c;
febio_spec.Material.material{1}.k1=k1;
febio_spec.Material.material{1}.k2=k2;
% febio_spec.Material.material{1}.gamma=gamma;
febio_spec.Material.material{1}.kappa=kappa;
febio_spec.Material.material{1}.k=k;

nElem = size(febio_spec.Mesh.Elements{1}.elem.VAL,1);

% Set orientation data and deformation magnitude 
for i = 1:1:numel(LangerAngles)

    % Create material axis vectors 
    R = euler2DCM([0,-LangerAngles(i),0]); % Rotation matrix
    e1 = (R*[1 0 0]')'; % Rotated e1 vector 
    e1_dir = e1(ones(nElem,1),:); % e1 copied for all elements        
    e2 = [0 1 0];% (R*[0 1 0]')'; % Rotated e1 vector 
    e2_dir = e2(ones(nElem,1),:); % e1 copied for all elements
    
    %Update the MeshData section with new material axis vectors 
    for q=1:1:nElem
        febio_spec.MeshData.ElementData{1}.elem{q}.ATTR.lid=q;
        febio_spec.MeshData.ElementData{1}.elem{q}.a=e1_dir(q,:);
        febio_spec.MeshData.ElementData{1}.elem{q}.d=e2_dir(q,:);
    end

    % Set the current displacement magnitude
    febio_spec.Boundary.bc{4}.value.VAL=displacementMagnitudes(i);

    % Export to .feb file
    febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
    
    % Run the current model 
    [runFlags(i)]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

    if runFlags(i) == 1
        % Importing nodal displacements from a log file
        simData(i).dataStruct_disp=importFEBio_logfile(febioLogFileName_disp,0,1);
    
        % Importing element stress from a log file
        simData(i).dataStruct_stress=importFEBio_logfile(febioLogFileName_stress,0,1);
    
        % Importing element stretch from a log file
        simData(i).dataStruct_stretch=importFEBio_logfile(febioLogFileName_stretch,0,1);
    else        
        simData(i).dataStruct_disp=NaN;
        simData(i).dataStruct_stress=NaN;
        simData(i).dataStruct_stretch=NaN;        
    end
end

if all(runFlags==1)
    stretch_perp_sim = squeeze(mean(mean(simData(1).dataStruct_stretch.data,1),2));
    stress_perp_sim = squeeze(mean(mean(simData(1).dataStruct_stress.data,1),2));
    
    stretch_para_sim = squeeze(mean(mean(simData(2).dataStruct_stretch.data,1),2));
    stress_para_sim = squeeze(mean(mean(simData(2).dataStruct_stress.data,1),2));
    
    % Interpolating simulation data for the experimental stretch points
    stress_perp_sim_exp = interp1(stretch_perp_sim,stress_perp_sim,stretch_exp_perp_raw,'pchip','extrap');
    stress_para_sim_exp = interp1(stretch_para_sim,stress_para_sim,stretch_exp_para_raw,'pchip','extrap');
    
    % Update graph
    figure(hf);
    set(h(1),'XData',stretch_perp_sim,'YData',stress_perp_sim);
    set(h(2),'XData',stretch_para_sim,'YData',stress_para_sim);

    set(h(3),'XData',stretch_exp_perp_raw,'YData',stress_perp_sim_exp);
    set(h(4),'XData',stretch_exp_para_raw,'YData',stress_para_sim_exp);

    % for q=1:1:numel(stress_perp_sim_exp)
    %     plotV([stretch_exp_perp_raw(q) stress_perp_sim_exp(q); stretch_exp_perp_raw(q) stress_exp_perp_raw(q);],'y-','LineWidth',2);
    % end
    % 
    % for q=1:1:numel(stress_para_sim_exp)
    %     plotV([stretch_exp_para_raw(q) stress_para_sim_exp(q); stretch_exp_para_raw(q) stress_exp_para_raw(q);],'c-','LineWidth',2);
    % end

    drawnow; axis tight; drawnow; 

    % Compute error metrics    
    diff1 = (stress_exp_perp_raw(2:end)-stress_perp_sim_exp(2:end)); % Difference
    errorVal1 = diff1; % errorVal1 = diff1./stress_exp_perp_raw(2:end); % Relative 

    diff2 = (stress_exp_para_raw(2:end)-stress_para_sim_exp(2:end)); % Difference
    errorVal2 = diff2; % errorVal2 = diff2./stress_exp_para_raw(2:end); % Relative 

    errorVal = [errorVal1(:); errorVal2(:)]; % Appending error vectors together     
else
    errorVal = [NaN(size(stress_exp_perp_raw)); NaN(size(stress_exp_para_raw))];
end
    
SSE = sum(errorVal.^2);
disp(['SSE: ', sprintf('%6.8e',SSE)]);

if objectiveStruct.method == 1 
    errorVal = SSE; % use sum of squares for Nelder-Mead
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
