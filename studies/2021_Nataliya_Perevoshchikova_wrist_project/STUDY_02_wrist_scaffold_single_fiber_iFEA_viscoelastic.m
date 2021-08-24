%% STUDY_02_wrist_scaffold_single_fiber_iFEA_viscoelastic
% Below is a demonstration for:
% 
% * Inverse FEA based material parameter optimisation 
% * febio_spec version 3.0
% * febio, FEBio
% * hexahedral elements, hex8
% * static, solid
% * hyperelastic, Ogden
% * displacement logfile
% * stress logfile


clear; close all; clc;

%%
% Plot settings
fontSize=40;
markerSize=25;
lineWidth=3; 

%% Control parameters
% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'2021_Nataliya_Perevoshchikova_wrist_project','data','temp');
loadPathMatProp=fullfile(defaultFolder,'2021_Nataliya_Perevoshchikova_wrist_project','data','mat_prop','PCL');

%%
% Building a quadrilateral circular mesh
pointSpacingHeight=0.1;
d=0.35;
h=10;
xSpacing=0.5;
ySpacing=0.4;
nCopies_x=8;
nCopies_y=8;

r=d/2;
ne=2; %Elements in radius
f=0.6; %Fraction (with respect to outer radius) where central square appears
NumStrands=4*7;
%Create the mesh
[Fq,Vq]=discQuadMesh(ne,r,f);

%%
% Lofting to a fiber
Vq(:,3)=0;
Vq1=Vq;
Vq1(:,3)=Vq1(:,3)-h/2;
Vq2=Vq;
Vq2(:,3)=Vq2(:,3)+h/2;

numStepsSweep=ceil(h/pointSpacingHeight);
xc=zeros(numStepsSweep,1);
yc=zeros(numStepsSweep,1);
zc=linspace(-h/2,h/2,numStepsSweep)';
Vc=[xc yc zc];

[~,~,~,S]=sweepLoft(Vq1,Vq2,[0 0 1],[0 0 1],Vc,numStepsSweep,0,0);

%Compose hexahedral elements
X=S.X'; Y=S.Y'; Z=S.Z'; %Coordinate matrices
V=[X(:) Y(:) Z(:)]; %Create node list

I=size(Vq,1)*((1:1:numStepsSweep)-1);
I=I(ones(size(Fq,1),1),:);
I=I(:);

FQ=repmat(Fq,numStepsSweep,1)+I(:,ones(size(Fq,2),1));
E=[FQ(1:end-size(Fq,1),:) FQ(size(Fq,1)+1:end,:)]; %The hexahedral elements
F_start=FQ(1:size(Fq,1),:);
F_end=FQ(end-size(Fq,1)+1:end,:);

[E,V,ind1,ind2]=mergeVertices(E,V); %Merge nodes (start and end are not shared yet) 

[F]=element2patch(E); %Element faces
F_start=ind2(F_start); %Start faces for BC's
F_end=ind2(F_end); %End faces for BC's
ind=tesBoundary(F,V); %Indices of boundary faces
Fb=F(ind,:); %Boundary faces (for plotting)


%%
% Visualizing mesh

cFigure; hold on;
title('Hexahedral mesh one fiber','fontSize',fontSize); 
gpatch(Fb,V,'gw','k',0.5);
gpatch(F_start,V,'rw','k',1);
gpatch(F_end,V,'bw','k',1);
axisGeom; 
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;


%%
% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress


C1=57.21894783089972;
C2=57.23071414361289;
K=11444.966197451262;
M1=2.00012552045284;
M2=-1.9997730268725347;
%Initial Viscoelastic parameters
g1_ini=1.25;%Viscoelastic QLV proportional coefficient
t1_ini=25;%Viscoelastic QLV time coefficient
P=[t1_ini g1_ini C1 M1 C2 M2];
displacementMagnitude=0.6030;

% Experimental viscoelastic properties
matname = fullfile(loadPathMatProp,'Viscoelastic_properties.mat');
load(matname);

cFigure; 
title('Viscoelastic region','fontSize',fontSize); 
hold on;
xlabel('Time, s','FontSize',fontSize);
ylabel('Load rate','FontSize',fontSize);
plot(V_time,lc,'k.-','LineWidth',3);
set(gca,'FontSize',fontSize);
camlight headlight;

[M,I]=max(lc);
time_load=V_time(1:I);
time_unload=V_time(I+1:end);
lc_load=lc(1:I);
lc_unload=lc(I+1:end);

% FEA control settings
max_refs=50; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=10; %Maximum number of retires
k_factor=100;
NumStrands=4*7;

t_load=max(time_load); %Time from start to max load
t_step_ini1=t_load/50; %Initial desired step size
numTimeSteps1=round(t_load/t_step_ini1); %Number of time steps desired
t_step1=t_load/numTimeSteps1; %Step size
dtmin1=t_step1/100; %Smallest allowed step size
dtmax1=t_step1; %Largest allowed step size

t_unload=V_time(end)-max(time_load);
t_step_ini2=t_step_ini1; %Initial desired step size
numTimeSteps2=round(t_unload/t_step_ini2); %Number of time steps desired
t_step2=t_unload/numTimeSteps2; %Step size
dtmin2=t_step2/100; %Smallest allowed step size
dtmax2=t_step2; %Largest allowed step size


cFigure; 
title('Viscoelastic region','fontSize',fontSize); 
hold on;
xlabel('Time, s','FontSize',fontSize);
ylabel('Load rate','FontSize',fontSize);
plot(time_load,lc_load,'r.-','LineWidth',3);
plot(time_unload,lc_unload,'k.-','LineWidth',3);
set(gca,'FontSize',fontSize);
camlight headlight;


%% DEFINE BC's
%Supported nodes
bcRigidList=unique(F_start(:));
 
%Prescribed force nodes
bcPrescribeList=unique(F_end(:));
bcPrescribeMagnitudes=displacementMagnitude(ones(1,numel(bcPrescribeList)),:);
 
%Visualize BC's
cFigure; hold on;
title('Boundary conditions','FontSize',fontSize);
gpatch(F,V,'y','k',1);

hl(1)=plotV(V(bcRigidList,:),'k.','MarkerSize',markerSize);
hl(2)=plotV(V(bcPrescribeList,:),'r.','MarkerSize',markerSize);

legend(hl,{'BC full support','BC prescribed Z displacement'})
axisGeom;
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow; 

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='3.0'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

analysisType='STATIC';
%Control sections for each step
febio_spec.Step.step{1}.Control=febio_spec.Control; %Copy from template
febio_spec.Step.step{1}.ATTR.id=1;
febio_spec.Step.step{1}.Control.analysis=analysisType;
febio_spec.Step.step{1}.Control.time_steps=numTimeSteps1;
febio_spec.Step.step{1}.Control.step_size=t_step1;
febio_spec.Step.step{1}.Control.solver.max_refs=max_refs;
febio_spec.Step.step{1}.Control.solver.max_ups=max_ups;
febio_spec.Step.step{1}.Control.time_stepper.dtmin=dtmin1;
febio_spec.Step.step{1}.Control.time_stepper.dtmax=dtmax1; 
febio_spec.Step.step{1}.Control.time_stepper.max_retries=max_retries;
febio_spec.Step.step{1}.Control.time_stepper.opt_iter=opt_iter;

febio_spec.Step.step{2}.Control=febio_spec.Control; %Copy from template
febio_spec.Step.step{2}.ATTR.id=2;
febio_spec.Step.step{2}.Control.analysis=analysisType;
febio_spec.Step.step{2}.Control.time_steps=numTimeSteps2;
febio_spec.Step.step{2}.Control.step_size=t_step2;
febio_spec.Step.step{2}.Control.solver.max_refs=max_refs;
febio_spec.Step.step{2}.Control.solver.max_ups=max_ups;
febio_spec.Step.step{2}.Control.time_stepper.dtmin=dtmin2;
febio_spec.Step.step{2}.Control.time_stepper.dtmax=dtmax2; 
febio_spec.Step.step{2}.Control.time_stepper.max_retries=max_retries;
febio_spec.Step.step{2}.Control.time_stepper.opt_iter=opt_iter;


%Remove control field (part of template) since step specific control sections are used
febio_spec=rmfield(febio_spec,'Control'); 

materialName1='Material1';
febio_spec.Material.material{1}.ATTR.name=materialName1;

%Viscoelastic part
febio_spec.Material.material{1}.ATTR.type='uncoupled viscoelastic';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.g1=g1_ini;
febio_spec.Material.material{1}.t1=t1_ini;

%Elastic part
febio_spec.Material.material{1}.elastic{1}.ATTR.type='Ogden';
febio_spec.Material.material{1}.elastic{1}.c1=C1;
febio_spec.Material.material{1}.elastic{1}.m1=M1;
febio_spec.Material.material{1}.elastic{1}.c2=C2;
febio_spec.Material.material{1}.elastic{1}.m2=M2;
febio_spec.Material.material{1}.elastic{1}.k=K;

    
%Geometry section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% % -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of the element set
febio_spec.Mesh.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type of this set
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E;
% 
 
%MeshDomains section
febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

% % -> NodeSets
% -> NodeSets
nodeSetName1='bcSupportList';
nodeSetName2='bcPrescribeList';

febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.node.ATTR.id=bcRigidList(:);
% 
febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.node.ATTR.id=bcPrescribeList(:);
% 
% %Boundary condition section 
% % -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.type='fix';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.dofs='x';

febio_spec.Boundary.bc{2}.ATTR.type='fix';
febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{2}.dofs='y';

febio_spec.Boundary.bc{3}.ATTR.type='fix';
febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{3}.dofs='z';

febio_spec.Boundary.bc{4}.ATTR.type='fix';
febio_spec.Boundary.bc{4}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{4}.dofs='x';

febio_spec.Boundary.bc{5}.ATTR.type='fix';
febio_spec.Boundary.bc{5}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{5}.dofs='y';

febio_spec.Boundary.bc{6}.ATTR.type='prescribe';
febio_spec.Boundary.bc{6}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{6}.dof='z';
febio_spec.Boundary.bc{6}.scale.ATTR.lc=1;
febio_spec.Boundary.bc{6}.scale.VAL=displacementMagnitude;
febio_spec.Boundary.bc{6}.relative=0;

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0;time_load lc_load;time_unload lc_unload];

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

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='sx';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.element_data{1}.VAL=1:size(E,1);

%% Quick viewing of the FEBio input file structure
% The |febView| function can be used to view the xml structure in a MATLAB
% figure window. 

%%
% febView(febio_spec); %Viewing the febio file|

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
febioAnalysis.runMode='internal';%'external';


[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!
  

%% Import FEBio results 

if runFlag==1 %i.e. a succesful run
    
    % Importing nodal displacements from a log file
    [~, N_disp_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp));
    
    N_disp_mat=N_disp_mat(:,2:end,:);
    sizImport=size(N_disp_mat);
    sizImport(3)=sizImport(3)+1;
    N_disp_mat_n=zeros(sizImport);
    N_disp_mat_n(:,:,2:end)=N_disp_mat;
    N_disp_mat=N_disp_mat_n;
    DN=N_disp_mat(:,:,end);
    DN_magnitude=sqrt(sum(DN(:,1).^2,2)); %X
    V_def=V+DN;
    [CF]=vertexToFaceMeasure(F,DN_magnitude);

    %Importaing force from log file
    [time_mat, Force_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_force));
    time_mat=[0; time_mat(:)]; %Time
    Force_mat=Force_mat(:,2:end,:);
    sizImport=size(Force_mat);
    sizImport(3)=sizImport(3)+1;
    Force_mat_n=zeros(sizImport);
    Force_mat_n(:,:,2:end)=Force_mat;
    Force_mat=Force_mat_n;
    time_sim=time_mat;
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp=gpatch(F,V_def,CF,'k',1); %Add graphics object to animate
    gpatch(F,V,0.5*ones(1,3),'none',0.25); %A static graphics object
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]);    
    axis([min(V_def(:,1)) max(V_def(:,1)) min(V_def(:,2)) max(V_def(:,2)) min(V_def(:,3)) max(V_def(:,3))]); %Set axis limits statically
    view(130,25); %Set view direction
    camlight headlight;        
        
    % Set up animation features
    animStruct.Time=time_mat; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        DN=N_disp_mat(:,:,qt); %Current displacement
        DN_magnitude=sqrt(sum(DN.^2,2)); %Current displacement magnitude
        V_def=V+DN; %Current nodal coordinates
        [CF]=vertexToFaceMeasure(F,DN_magnitude); %Current color data to use
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,CF}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;
    

    %% 
    % Calculate the simulated applied uniaxial stretch
    DZ_set=N_disp_mat(bcPrescribeList,3,:); %Z displacements of the prescribed set
    DZ_set=mean(DZ_set,1); %Calculate mean Z displacements across nodes
    DZ_set=squeeze(DZ_set);
    
    %Calculated the simulated force
    FZ_set=Force_mat(bcPrescribeList,3,:);
    FZ_set=sum(FZ_set,1); %Calculate mean X force across nodes
    FZ_set=squeeze(FZ_set);
    
    disp_sim=DZ_set;
    load_sim=FZ_set*NumStrands;
    
        %%    
    % Visualize curves
 
    hf=cFigure;
    hold on;    
    title('Cyclic load-displacement curve','FontSize',fontSize);
    xlabel('Time','FontSize',fontSize); ylabel('Load','FontSize',fontSize); 

    

    hl(1)=plot(time_mat,load_sim,'r.-','lineWidth',lineWidth*1.1,'markerSize',markerSize*1.3);
    hl(2)=plot(V_time,V_load,'k-','lineWidth',lineWidth);
    
    legend(hl,{'Simulation','Experiment'},'Location','northwest');
    

    view(2); axis tight;  grid on; axis square; box on; 
    set(gca,'FontSize',fontSize);
    drawnow;

end

%% Create structures for optimization 
% Material structure
mat_struct.id=1; %Material id


mat_struct.par_names={'t1','g1','c1','m1','c2','m2','k'}; %Parameter names
mat_struct.par_values={t1_ini g1_ini C1 M1 C2 M2 K}; %Parameter values

objectiveStruct.parNormFactors=abs(P); %This will normalize the parameters to ones(size(P))
objectiveStruct.Pb_struct.xx_c=P; %Parameter constraining centre
objectiveStruct.Pb_struct.xxlim=[[P(1)/100 0.01 P(3)/100 -30 P(5)/100 -30]' [P(1)*1e3 50 P(3)*100 30 P(5)*100 30]']; %Parameter bounds
 

febioAnalysis.time_on=0; 
febioAnalysis.time_log_on=0; 

%What should be known to the objective function:
objectiveStruct.h=hl(1);
objectiveStruct.bcPrescribeList=bcPrescribeList;
objectiveStruct.time_exp=V_time;
objectiveStruct.load_exp=V_load;
objectiveStruct.febioAnalysis=febioAnalysis;
objectiveStruct.febio_spec=febio_spec;
objectiveStruct.febioFebFileName=febioFebFileName;
objectiveStruct.mat_struct=mat_struct;
objectiveStruct.k_factor=k_factor;
objectiveStruct.NumStrands=NumStrands;


%Optimisation settings
maxNumberIterations=100; %Maximum number of optimization iterations
maxNumberFunctionEvaluations=maxNumberIterations*10; %Maximum number of function evaluations, N.B. multiple evaluations are used per iteration
functionTolerance=1e-25; %Tolerance on objective function value
parameterTolerance=0.001; %Tolerance on parameter variation
displayTypeIterations='iter';

objectiveStruct.method=2; 

%File names of output files
output_names.load=fullfile(savePath,febioLogFileName_force);

objectiveStruct.run_output_names=output_names;

%% start optimization

Pn=P./objectiveStruct.parNormFactors;

switch objectiveStruct.method
    case 1 %fminsearch and Nelder-Mead
        OPT_options=optimset('fminsearch'); % 'Nelder-Mead simplex direct search'
        OPT_options = optimset(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
            'MaxIter',maxNumberIterations,...
            'TolFun',functionTolerance,...
            'TolX',parameterTolerance,...
            'Display',displayTypeIterations,...
            'FinDiffRelStep',1e-2,...
            'DiffMaxChange',0.5);
        [Pn_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= fminsearch(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct,savePath),Pn,OPT_options);
    case 2 %lsqnonlin and Levenberg-Marquardt
        OPT_options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
        OPT_options = optimoptions(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
            'MaxIter',maxNumberIterations,...
            'TolFun',functionTolerance,...
            'TolX',parameterTolerance,...
            'Display',displayTypeIterations,...
            'FinDiffRelStep',1e-2,...
            'DiffMaxChange',0.5);
        [Pn_opt,OPT_out.resnorm,OPT_out.residual]= lsqnonlin(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct,savePath),Pn,[],[],OPT_options);    
end

%%
[Fopt,OPT_stats_out]=objectiveFunctionIFEA(Pn_opt,objectiveStruct,savePath);

%% Unnormalize and constrain parameters

P_opt=Pn_opt.*objectiveStruct.parNormFactors; %Scale back, undo normalization

%Constraining parameters
for q=1:1:numel(P)
    [P(q)]=boxconstrain(P(q),objectiveStruct.Pb_struct.xxlim(q,1),objectiveStruct.Pb_struct.xxlim(q,2),objectiveStruct.Pb_struct.xx_c(q));    
end

disp_text=sprintf('%6.16e,',P_opt); disp_text=disp_text(1:end-1);
disp(['P_opt=',disp_text]);


%%
hf1=cFigure;
%title('displacement load curves optimised','FontSize',fontSize);
xlabel('Time (s)','FontSize',fontSize); ylabel('Force (N)','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
   
Hn(1)=plot(OPT_stats_out.time_sim,OPT_stats_out.load_sim,'r.-','lineWidth',lineWidth*1.1,'markerSize',markerSize*1.3);
%Hn(1)=plot(OPT_stats_out.time_sim,OPT_stats_out.load_sim,'r-','lineWidth',lineWidth,'markerSize',markerSize);
%Hn(2)=plot(V_time,V_load,'g-','lineWidth',lineWidth);
Hn(2)=plot(V_time,V_load,'k-','lineWidth',lineWidth);

legend(Hn,{'Simulation','Experiment'},'Location','northwest');
view(2); axis tight;  grid on;
set(gca,'FontSize',fontSize);
drawnow;


function [Fopt,OPT_stats_out]=objectiveFunctionIFEA(Pn,objectiveStruct,savePath)

%%

febioFebFileName=objectiveStruct.febioFebFileName; 
febio_spec=objectiveStruct.febio_spec; 

%% Unnormalize and constrain parameters

P=Pn.*objectiveStruct.parNormFactors; %Scale back, undo normalization
P_in=P; %Proposed P


%Constraining parameters
for q=1:1:numel(P)
    [P(q)]=boxconstrain(P(q),objectiveStruct.Pb_struct.xxlim(q,1),objectiveStruct.Pb_struct.xxlim(q,2),objectiveStruct.Pb_struct.xx_c(q));    
end

%% Setting material parameters

%Acces material parameters

mat_struct=objectiveStruct.mat_struct;
mat_struct.par_values={P(1) P(2) P(3) P(4) P(5) P(6) (P(3)+P(5))*objectiveStruct.k_factor};
   

disp('SETTING MATERIAL PARAMETERS...');
disp(['Proposed (norm.): ',sprintf(repmat('%6.6e ',[1,numel(Pn)]),Pn)]);
disp(['Proposed        : ',sprintf(repmat('%6.6e ',[1,numel(P_in)]),P_in)]);
disp(['Set (constr.)   : ',sprintf(repmat('%6.6e ',[1,numel(P)]),P)]);

%Assign material parameters
matId=mat_struct.id;
for q=1:1:2
    parNameNow=mat_struct.par_names{q};
    parValuesNow=mat_struct.par_values{q};    
    febio_spec.Material.material{matId}.(parNameNow)=mat2strIntDouble(parValuesNow);
end

for q=3:1:numel(mat_struct.par_names)
    parNameNow=mat_struct.par_names{q};
    parValuesNow=mat_struct.par_values{q};    
    febio_spec.Material.material{matId}.elastic{1}.(parNameNow)=mat2strIntDouble(parValuesNow);
end

%febView(febio_spec); %Viewing the febio file|

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
disp('Done')

%% START FEBio

[runFlag]=runMonitorFEBio(objectiveStruct.febioAnalysis);

pause(0.1); 

bcPrescribeList=objectiveStruct.bcPrescribeList;
time_exp=objectiveStruct.time_exp;
load_exp=objectiveStruct.load_exp;

if runFlag==1    

    %Importing force from log file
    [time_mat, Force_mat,~]=importFEBio_logfile(objectiveStruct.run_output_names.load);
    time_mat=[0; time_mat(:)]; %Time
    Force_mat=Force_mat(:,2:end,:); %Remove index
    sizImport=size(Force_mat);
    sizImport(3)=sizImport(3)+1;
    Force_mat_n=zeros(sizImport);
    Force_mat_n(:,:,2:end)=Force_mat; 
    Force_mat=Force_mat_n; %With added zeros for t=0
        
    %Calculated the simulated force
    FZ_set=Force_mat(bcPrescribeList,3,:); 
    FZ_set=sum(FZ_set,1); %Calculate mean X force across nodes
    FZ_set=squeeze(FZ_set);

    time_sim=time_mat;
    load_sim=FZ_set*objectiveStruct.NumStrands; %Multiply force by number of strands
    
    if ~isempty(objectiveStruct.h)
        objectiveStruct.h.XData=time_sim;
        objectiveStruct.h.YData=load_sim;
        drawnow;
    end

    load_sim_exp = interp1(time_sim,load_sim,time_exp,'pchip');

    %Derive Fopt
    loadDev=load_exp-load_sim_exp;
   
    switch objectiveStruct.method
        case 1
            Fopt=sum((loadDev).^2)/length(loadDev); %Sum of squared differences
        case 2
            Fopt=sum((loadDev).^2)/length(loadDev);%loadDev(:); %Squared differences
    end

    OPT_stats_out.load_sim=load_sim;
    OPT_stats_out.time_sim=time_sim;
    OPT_stats_out.loadDev=loadDev;
    OPT_stats_out.Fopt=Fopt;
    
    
else %Output NaN
    switch objectiveStruct.method
        case 1
            Fopt=NaN; 
        case 2
            Fopt=NaN(size(load_exp)); %Squared differences
    end
    OPT_stats_out=[];
end

end