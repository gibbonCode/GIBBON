%% DEMO_create_run_import_FEBIO_spheres
% Below is a demonstration for:
% 
% * The use of TETgen for meshing based on surface geometry
% * The specification of boundary conditions for FEBio
% * The exporting of .feb files
% * Running an FEBio job with MATLAB
% * Importing FEBio results into MATLAB

%%

clear; close all; clc; 

%%
% Plot settings
fontSize=15;
faceAlpha1=0.5;
faceAlpha2=0.5;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;

% path names
filePath=mfilename('fullpath');
savePath=fullfile(fileparts(filePath),'data','temp');

%% Defining the surface models
% The model will consists of two spheres one contained within the other
% defining two material regions. A stiff core and a soft outer later.

%%
% Control parameters for surface models
r1=2.75; %Outer sphere radius
numRefine1=3; %Number of refinement steps from icosahedron
faceBoundMarker1=2; %Face marker for outer sphere

r2=2; %Inner sphere radius
numRefine2=3; %Number of refinement steps from icosahedron
faceBoundMarker2=3; %Face marker for inner sphere

nSteps=20;
max_refs=25;
max_ups=15;

contactType=0; 
    
%%
% Building the spheres using |geoSphere| function

[F1,V1,~]=geoSphere(numRefine1,r1);
% F1=fliplr(F1);
[F2,V2,~]=geoSphere(numRefine2,r2);
% F2=fliplr(F2);

% Merging the model geometries into a single set
V=[V1;V2]; %Joining nodes
F=[F1;F2+size(V1,1)]; %Joining faces
faceBoundaryMarker=[faceBoundMarker1*ones(size(F1,1),1); faceBoundMarker2*ones(size(F2,1),1)]; %Create boundary markers for faces

%%
% Plotting surface models
hf=cFigure;
title('Surface models','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',faceBoundaryMarker,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
[hp]=patchNormPlot(F,V,0.25);

colormap(autumn(2));
colorbar;
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;


%% CREATING A SOLID TETRAHEDRAL MESH USING TETGEN

%%
% First region points need to be defined. These represent a list of
% arbitrary coordinates for points inside the regions. 1 point per region
% is specified.
% For the example here the points are easily specified. Sometimes a
% raytracing algorythm or the use of the |triSurf2Im| function is required
% to find interior points.

V_regions=[0 0 (r1+r2)/2;0 0 0;]; % Define region points

%%
% Next holes are defined. These are similar to regions. However holes, as
% the name suggests, are regions that a not meshed and are left empty.
% This model does not contain holes so the list is empty

V_holes=[]; %Define hole points

%%
% For each region the mesh density parameter can be specified
[v]=tetVolMeanEst(F,V); %Estimate volume of ideal tetrahedron
regionA=[v v]; % Regional mesh parameters

%%
% CREATING THE SMESH STRUCTURE.
% TetGen can mesh geometries from various mesh file formats. For the GIBBON
% toolbox .smesh files have been implemented. Below a structure is created
% that fully defines such as smesh file and the meshing settings for
% TetGen.

stringOpt='-pq1.2AaYQ';
modelNameEnd='tempModel';
modelName=fullfile(savePath,modelNameEnd);

inputStruct.stringOpt=stringOpt;
inputStruct.Faces=F;
inputStruct.Nodes=V;
inputStruct.holePoints=V_holes;
inputStruct.faceBoundaryMarker=faceBoundaryMarker; %Face boundary markers
inputStruct.regionPoints=V_regions; %region points
inputStruct.regionA=regionA;
inputStruct.minRegionMarker=2; %Minimum region marker
inputStruct.modelName=modelName;

%%
% Mesh model using tetrahedral elements using tetGen (see:
% <http://wias-berlin.de/software/tetgen/>)

[meshOutput]=runTetGen(inputStruct); %Run tetGen

%%
% Accessing the model element and patch data
FT=meshOutput.faces;
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
VT=meshOutput.nodes;
C=meshOutput.faceMaterialID;
E=meshOutput.elements;
elementMaterialIndices=meshOutput.elementMaterialID;

%%
% Plotting the meshed geometry

hf1=cFigure;
subplot(1,3,1);
title('Solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',FT,'Vertices',VT,'FaceColor','flat','CData',C,'lineWidth',edgeWidth,'edgeColor',edgeColor);
% [hp]=patchNormPlot(FT,VT,0.25);
view(3); axis tight;  axis equal;  grid on;
colormap(autumn);
% camlight headlight;
set(gca,'FontSize',fontSize);

subplot(1,3,2);
title('Model boundaries','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',Fb,'Vertices',VT,'FaceColor','flat','CData',Cb,'lineWidth',edgeWidth,'edgeColor',edgeColor,'FaceAlpha',faceAlpha1);
view(3); axis tight;  axis equal;  grid on;
colormap(autumn);
set(gca,'FontSize',fontSize);
drawnow;

subplot(1,3,3);
%Selecting half of the model to see interior
Y=VT(:,2); YE=mean(Y(E),2);
L=YE>mean(Y);
[Fs,Cs]=element2patch(E(L,:),C(L));

title('Cut view of solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',Fs,'Vertices',VT,'FaceColor','flat','CData',Cs,'lineWidth',edgeWidth,'edgeColor',edgeColor);

% [hp]=patchNormPlot(Fs,VT,0.25);

view(3); axis tight;  axis equal;  grid on;
colormap(autumn);
% camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;

%% Splitting up material regions

Em1=E(elementMaterialIndices==-2,:);
ind_m=unique(Em1(:));
Vm1=VT(ind_m,:);
indFix=1:size(VT,1);
indFix(ind_m)=1:numel(ind_m);
Em1=indFix(Em1);

F1=Fb(Cb==3,:);
F1=indFix(F1);
F1=fliplr(F1);

Em2=E(elementMaterialIndices==-3,:);
ind_m=unique(Em2(:));
Vm2=VT(ind_m,:);
indFix=1:size(VT,1);
indFix(ind_m)=1:numel(ind_m);
Em2=indFix(Em2);

F2=Fb(Cb==3,:);
F2=indFix(F2);

%% Joining node sets

VT=[Vm1;Vm2];
Em2=Em2+size(Vm1,1);
F2=F2+size(Vm1,1);

%%

[Fm1]=element2patch(Em1);
[Fm2]=element2patch(Em2);

cFigure;
subplot(1,2,1); hold on;
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;

hps=patch('Faces',Fm1,'Vertices',VT,'FaceColor','b','lineWidth',edgeWidth,'edgeColor',edgeColor);

view(3); axis tight;  axis equal;  grid on;
colormap(autumn);
camlight headlight;
set(gca,'FontSize',fontSize);

subplot(1,2,2); hold on;
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;

hps=patch('Faces',Fm2,'Vertices',VT,'FaceColor','r','lineWidth',edgeWidth,'edgeColor',edgeColor);

view(3); axis tight;  axis equal;  grid on;
colormap(autumn);
camlight headlight;
set(gca,'FontSize',fontSize);

%%
cFigure; hold on;
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;

hps=patch('Faces',F1,'Vertices',VT,'FaceColor','b','lineWidth',edgeWidth,'edgeColor',edgeColor,'FaceAlpha',0.25);
[hp]=patchNormPlot(F1,VT,0.25);
hps=patch('Faces',F2,'Vertices',VT,'FaceColor','r','lineWidth',edgeWidth,'edgeColor',edgeColor,'FaceAlpha',0.25);
[hp]=patchNormPlot(F2,VT,0.25);

view(3); axis tight;  axis equal;  grid on;
colormap(autumn);
% camlight headlight;
set(gca,'FontSize',fontSize);

drawnow; 

%% DEFINE PRESCRIBED DISPLACEMENTS
% For this example the outer sphere nodes are subjected to a prescribed displacement

%Get outer surface (numbering may have altered due to tetgen behaviour so
%redefined here)
indOuter=unique(F1(:));
V1=VT(indOuter,:);

[N,Vn,Nv]=patchNormal(F1,VT);

ampDef=0.5; 

V1_def=V1+Nv(indOuter,:)*ampDef;


% Define boundary displacement values
bcPrescribedMagnitudes=(V1_def-V1);

% Define indices (node numbers) for the prescribed displacement
bcIndicesPrescribed=indOuter;

%%
% Plotting deformed outer surface
C_outer=sqrt(sum(bcPrescribedMagnitudes.^2,2)); %Color towards displacement magnitude
CV=zeros(size(VT,1),1);
CV(indOuter)=C_outer;
[CF]=vertexToFaceMeasure(F1,CV);

VT_def=VT;
VT_def(indOuter,:)=V1_def;

hf=cFigure;
title('The deformed outer surface','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F1,'Vertices',VT_def,'FaceColor','flat','CData',CF,'FaceAlpha',1);
colormap jet; colorbar;
% camlight headlight;
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
FEB_struct.Geometry.Nodes=VT;
FEB_struct.Geometry.Elements={Em1,Em2}; %The element sets
FEB_struct.Geometry.ElementType={'tet4','tet4'}; %The element types
FEB_struct.Geometry.ElementMat={ones(1,size(Em1,1)) 2*ones(1,size(Em2,1))};
FEB_struct.Geometry.ElementsPartName={'Sphere_1','Sphere_2'};

%Defining surfaces
FEB_struct.Geometry.Surface{1}.Set=F1;
FEB_struct.Geometry.Surface{1}.Type='tri3';
FEB_struct.Geometry.Surface{1}.Name='Contact_master';

FEB_struct.Geometry.Surface{2}.Set=F2;
FEB_struct.Geometry.Surface{2}.Type='tri3';
FEB_struct.Geometry.Surface{2}.Name='Contact_slave';

%Adding contact information
FEB_struct.Contact{1}.Surface{1}.SetName=FEB_struct.Geometry.Surface{1}.Name;
FEB_struct.Contact{1}.Surface{1}.Type='master';
FEB_struct.Contact{1}.Surface{2}.SetName=FEB_struct.Geometry.Surface{2}.Name;
FEB_struct.Contact{1}.Surface{2}.Type='slave';
switch contactType
    case 0 %Tied
        FEB_struct.Contact{1}.Type='tied';
        FEB_struct.Contact{1}.Properties={'penalty','tolerance'};
        FEB_struct.Contact{1}.Values={50,0.1};
    case 1 %STICKY
        FEB_struct.Contact{1}.Type='sticky';
        FEB_struct.Contact{1}.Properties={'laugon','tolerance','penalty','minaug','maxaug'};
        FEB_struct.Contact{1}.Values={0,0.1,100,0,10};
    case 2 %SLIDING facet-facet
        FEB_struct.Contact{1}.Type='facet-to-facet sliding';
        FEB_struct.Contact{1}.Properties={'penalty','auto_penalty','two_pass',...
            'laugon','tolerance',...
            'gaptol','minaug','maxaug',...
            'seg_up',...
            'search_tol'};
        FEB_struct.Contact{1}.Values={50,1,0,...
            0,0.1,...
            0,0,10,...
            0,...
            0.01};
    case 3 %SLIDING sliding_with_gaps
        FEB_struct.Contact{1}.Type='sliding_with_gaps';
        FEB_struct.Contact{1}.Properties={'penalty','auto_penalty','two_pass',...
            'laugon','tolerance',...
            'gaptol','minaug','maxaug',...
            'fric_coeff','fric_penalty',...
            'seg_up',...
            'search_tol'};
        FEB_struct.Contact{1}.Values={50,1,0,...
            0,0.1,...
            0,0,10,...
            0,1,...
            0,...
            0.01};
    case 4 %sliding-tension-compression
        FEB_struct.Contact{1}.Type='sliding-tension-compression';
        FEB_struct.Contact{1}.Properties={'penalty','auto_penalty','two_pass',...
            'laugon','tolerance',...
            'gaptol','minaug','maxaug',...
            'fric_coeff','fric_penalty',...
            'seg_up',...
            'search_tol'};
        FEB_struct.Contact{1}.Values={50,1,0,...
            0,0.1,...
            0,0,10,...
            0,1,...
            0,...
            0.01};
end

% DEFINING MATERIALS

c1=1e-3;
k=c1*1e3;

%Material 1 Homes-Mow compressible multigeneration
c1_g=[c1/1000 c1/10];
k_g=1e4*c1_g;
FEB_struct.Materials{1}.Type='multigeneration';
FEB_struct.Materials{1}.Name='sphere_skin_mat';
FEB_struct.Materials{1}.Generation{1}.Solid{1}.Type='coupled Mooney-Rivlin';
FEB_struct.Materials{1}.Generation{1}.Solid{1}.Properties={'c1','c2','k'};
FEB_struct.Materials{1}.Generation{1}.Solid{1}.Values={c1_g(1),0,k_g(1)};
FEB_struct.Materials{1}.Generation{1}.Properties={'start_time'};
FEB_struct.Materials{1}.Generation{1}.Values={0};
FEB_struct.Materials{1}.Generation{2}.Solid{1}.Type='coupled Mooney-Rivlin';
FEB_struct.Materials{1}.Generation{2}.Solid{1}.Properties={'c1','c2','k'};
FEB_struct.Materials{1}.Generation{2}.Solid{1}.Values={c1_g(2),0,k_g(2)};
FEB_struct.Materials{1}.Generation{2}.Properties={'start_time'};
FEB_struct.Materials{1}.Generation{2}.Values={1};

%Material 2 deformable block
FEB_struct.Materials{2}.Type='coupled Mooney-Rivlin';
FEB_struct.Materials{2}.Name='sphere_mat';
FEB_struct.Materials{2}.Properties={'c1','c2','k'};
FEB_struct.Materials{2}.Values={c1,0,k};

%Defining node sets
FEB_struct.Geometry.NodeSet{1}.Set=bcIndicesPrescribed;
FEB_struct.Geometry.NodeSet{1}.Name='set_1';

%Adding BC information
FEB_struct.Boundary.Prescribe{1}.Set=bcIndicesPrescribed;
FEB_struct.Boundary.Prescribe{1}.bc='x';
FEB_struct.Boundary.Prescribe{1}.lc=1;
FEB_struct.Boundary.Prescribe{1}.nodeScale=bcPrescribedMagnitudes(:,1);
FEB_struct.Boundary.Prescribe{2}.Set=bcIndicesPrescribed;
FEB_struct.Boundary.Prescribe{2}.bc='y';
FEB_struct.Boundary.Prescribe{2}.lc=1;
FEB_struct.Boundary.Prescribe{2}.nodeScale=bcPrescribedMagnitudes(:,2);
FEB_struct.Boundary.Prescribe{3}.Set=bcIndicesPrescribed;
FEB_struct.Boundary.Prescribe{3}.bc='z';
FEB_struct.Boundary.Prescribe{3}.lc=1;
FEB_struct.Boundary.Prescribe{3}.nodeScale=bcPrescribedMagnitudes(:,3);

%Adding output requests
FEB_struct.Output.VarTypes={'displacement','stress','relative volume','contact pressure','contact traction','contact force'};

%Specify log file output
run_node_output_name=[modelNameEnd,'_node_out.txt'];
FEB_struct.run_output_names={run_node_output_name};
FEB_struct.output_types={'node_data'};
FEB_struct.data_types={'ux;uy;uz'};

%Step specific control sections
FEB_struct.Step{1}.Control.AnalysisType='static';
FEB_struct.Step{1}.Control.Properties={'time_steps','step_size',...
    'max_refs','max_ups',...
    'dtol','etol','rtol','lstol'};

tStep=1/nSteps;
FEB_struct.Step{1}.Control.Values={nSteps,tStep,...
    max_refs,max_ups,...
    0.001,0.01,0,0.9};
FEB_struct.Step{1}.Control.TimeStepperProperties={'dtmin','dtmax','max_retries','opt_iter','aggressiveness'};
FEB_struct.Step{1}.Control.TimeStepperValues={tStep/100,tStep,5,round(max_refs/2), 1};

FEB_struct.Step{2}.Control=FEB_struct.Step{1}.Control;

%Load curves
FEB_struct.LoadData.LoadCurves.id=1;
FEB_struct.LoadData.LoadCurves.type={'linear'};
FEB_struct.LoadData.LoadCurves.loadPoints={[0 0;1 1; 2 -0.1;]};

%% SAVING .FEB FILE

FEB_struct.disp_opt=0; %Display waitbars
febStruct2febFile(FEB_struct);

%% RUNNING FEBIO JOB

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
    [FEBIO_savePath,~,~] = fileparts(FEB_struct.run_filename); 
    [~, N_disp_mat,~]=importFEBio_logfile(fullfile(FEBIO_savePath,FEB_struct.run_output_names{1})); %Nodal displacements
    
    DN=N_disp_mat(:,2:end,end); %Final nodal displacements
    
    %% CREATING NODE SET IN DEFORMED STATE
    VT_def=VT+DN;
    
    %%
    % Plotting the meshed geometry
    
    %Selecting half of the model to see interior
    Z=VT(:,3); ZE=mean(Z(E),2);
    L=ZE<mean(Z);
    [Fs,~]=element2patch(E(L,:),[]);
    
    Cs=sqrt(sum(DN.^2,2)); %Color towards displacement magnitude
    
    hf1=cFigure;
    title('Cut view of deformed model showing internal results','FontSize',fontSize);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
    
    hps=patch('Faces',Fs,'Vertices',VT_def,'FaceColor','flat','FaceVertexCData',Cs);
    
    view(3); axis tight;  axis equal;  grid on;
    colormap jet; colorbar; shading interp;
    % camlight headlight;
    set(gca,'FontSize',fontSize);
    drawnow;
    
    %% EXAMPLE FOR VISUALIZATION OF MODEL OUTER SURFACE ONLY
    % Visualizing the outer surface only is less memory intensive for large models
    
    %Get free faces
    TR = triangulation(E,VT_def); %"Triangulation" representation
    F_free = freeBoundary(TR); %Free boundary triangles i.e. outer surface
    ind_V_free =unique(F_free(:)); %Indices of nodes at free boundary
    
    %Compute an example distance metric for visualization
    D=minDist(VT_def(ind_V_free,:),VT(ind_V_free,:));
    
    %Disance metric is known for a list of points not suitable yet for colouring
    %faces
    C=zeros(size(VT,1),1); %Initialse vertex color list
    C(ind_V_free)=D; %Set color for point selection
    [CF]=vertexToFaceMeasure(F_free,C); %Convert vertex to face color measure
    
    hf1=cFigure;
    title('Outer surface only with distance metric','FontSize',fontSize);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
    
    hps=patch('Faces',F_free,'Vertices',VT_def,'FaceColor','flat','CData',CF);
    hps=patch('Faces',F_free,'Vertices',VT,'FaceColor',0.5.*ones(1,3),'FaceAlpha',0.5,'EdgeColor','none');
    
    view(3); axis tight;  axis equal;  grid on;
    colormap jet; colorbar;
    % camlight headlight;
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