%% DEMO_create_run_import_FEBIO_spheres_pressure_v1p2
% Below is a demonstration for: 
% 
% *  The use of TETgen for meshing based on surface geometry
% *  The specification of boundary conditions for FEBio including pressure
% loading
% *  The exporting of .feb files
% *  Running an FEBio job with MATLAB
% *  Importing FEBio results into MATLAB

%%

clear; close all; clc;

%%
% Plot settings
figColor='w'; figColorDef='white';
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
r1=2; %Outer sphere radius
numRefine1=2; %Number of refinement steps from icosahedron
faceBoundMarker1=2; %Face marker for outer sphere

r2=1.25; %Inner sphere radius
numRefine2=2; %Number of refinement steps from icosahedron
faceBoundMarker2=3; %Face marker for inner sphere

%%
% Building the spheres using |geoSphere| function

[F1,V1,~]=geoSphere(numRefine1,r1);
[F2,V2,~]=geoSphere(numRefine2,r2);

% Merging the model geometries into a single set
V=[V1;V2]; %Joining nodes
F=[F1;F2+size(V1,1)]; %Joining faces
faceBoundaryMarker=[faceBoundMarker1*ones(size(F1,1),1); faceBoundMarker2*ones(size(F2,1),1)]; %Create boundary markers for faces

%%
% Plotting surface models
hf=figuremax(figColor,figColorDef);
title('Surface models','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',faceBoundaryMarker,'FaceAlpha',faceAlpha1,'lineWidth',edgeWidth,'edgeColor',edgeColor);
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
% regionA=[0.005 0.005]; % Regional mesh parameters

[edgeLengths]=patchEdgeLengths(F1,V1);
edgeLengthsMean=mean(edgeLengths);
meanProposedVolume=edgeLengthsMean^3./(6*sqrt(2)); %For regular tetrahedron
region1_A=0.5.*meanProposedVolume;

[edgeLengths]=patchEdgeLengths(F2,V2);
edgeLengthsMean=mean(edgeLengths);
meanProposedVolume=edgeLengthsMean^3./(6*sqrt(2)); %For regular tetrahedron
region2_A=0.5.*meanProposedVolume;

regionA=[region1_A region2_A];

%%
% CREATING THE SMESH STRUCTURE. 
% TetGen can mesh geometries from various mesh file formats. For the GIBBON
% toolbox .smesh files have been implemented. Below a structure is created
% that fully defines such as smesh file and the meshing settings for
% TetGen. 

stringOpt='-pq1.2AaYQ';
modelName=fullfile(savePath,'tempModel');
smeshName=[modelName,'.smesh'];

smeshStruct.stringOpt=stringOpt;
smeshStruct.Faces=F;
smeshStruct.Nodes=V;
smeshStruct.holePoints=V_holes;
smeshStruct.faceBoundaryMarker=faceBoundaryMarker; %Face boundary markers
smeshStruct.regionPoints=V_regions; %region points
smeshStruct.regionA=regionA;
smeshStruct.minRegionMarker=2; %Minimum region marker
smeshStruct.smeshName=smeshName;

%% 
% Mesh model using tetrahedral elements using tetGen (see:
% <http://wias-berlin.de/software/tetgen/>)

[meshOutput]=runTetGenSmesh(smeshStruct); %Run tetGen

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

hf1=figuremax(figColor,figColorDef);
subplot(1,3,1);
title('Solid tetrahedral mesh model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
hps=patch('Faces',FT,'Vertices',VT,'FaceColor','flat','CData',C,'lineWidth',edgeWidth,'edgeColor',edgeColor);
view(3); axis tight;  axis equal;  grid on;
colormap(autumn);
camlight headlight;
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
view(3); axis tight;  axis equal;  grid on;
colormap(autumn);
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;

%% DEFINE FACES FOR PRESSURE
% For this example the outer sphere nodes are subjected to a pressure

%Get outer surface (numbering may have altered due to tetgen behaviour so
%redefined here)
F1=Fb(Cb==2,:); 
F1=fliplr(F1);

%% 

hf=figuremax(figColor,figColorDef);
title('The outer surface','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F1,'Vertices',VT,'FaceColor','flat','CData',Cb(Cb==2),'FaceAlpha',0.5);
[hp]=patchNormPlot(F1,VT,0.5);
colormap jet; colorbar; 
camlight headlight;
set(gca,'FontSize',fontSize);
view(3); axis tight;  axis equal;  grid on;

%% Define node numbers for rigid support

[minVal,minInd]=min(sum((VT-(ones(size(VT,1),1)*mean(VT,1))).^2,2)); %Find node closest to centre
boundaryConditionNodeList=minInd; %Node closest to centre is constrained from moving

%% CONSTRUCTING FEB MODEL

% Defining file names
FEB_struct.run_filename=[modelName,'.feb']; %FEB file name
FEB_struct.run_logname=[modelName,'.txt']; %FEBio log file name

febMatID=elementMaterialIndices;
febMatID(elementMaterialIndices==-2)=1;
febMatID(elementMaterialIndices==-3)=2;

%Creating FEB_struct
FEB_struct.Geometry.Nodes=VT;
FEB_struct.Geometry.Elements={E}; %The element sets
FEB_struct.Geometry.ElementType={'tet4'}; %The element types
FEB_struct.Geometry.ElementMat={febMatID};

% DEFINING MATERIALS
k_factor=1000;

%Material 1
c1=1e-3;
k=1;%c1*k_factor;
Mat1.type='Mooney-Rivlin';
Mat1.props={'c1','c2','k'};
Mat1.vals={c1,0,k};
Mat1.aniso_type='none';

%Material 2
c1=10e-3;
k=c1*k_factor;
Mat2.type='Mooney-Rivlin';
Mat2.props={'c1','c2','k'};
Mat2.vals={c1,0,k};
Mat2.aniso_type='none';
 
% %Material 1
% c1=1e-3;
% m1=12;
% cp=c1*k_factor;
% Mat1.type='Ogden unconstrained';
% Mat1.props={'c1','m1','cp'};
% Mat1.vals={c1,m1,cp};
% Mat1.aniso_type='none';
% 
% %Material 2
% c1=10e-3;
% m1=12;
% cp=c1*k_factor;
% Mat2.type='Ogden unconstrained';
% Mat2.props={'c1','m1','cp'};
% Mat2.vals={c1,m1,cp};
% Mat2.aniso_type='none';

FEB_struct.Materials{1}=Mat1;
FEB_struct.Materials{2}=Mat2;

%Adding BC information
FEB_struct.Boundary.FixList={boundaryConditionNodeList};
FEB_struct.Boundary.FixType={'xyz'};

%Adding load information
FEB_struct.Loads.Pressure.Elements={F1};
FEB_struct.Loads.Pressure.Scale={0.5*ones(size(F1,1),1)};
FEB_struct.Loads.Pressure.LoadCurveIds=[1];

%Adding output requests
FEB_struct.Output.VarTypes={'displacement','stress','relative volume','shell thickness'};

%Specify log file output
run_node_output_name=[FEB_struct.run_filename(1:end-4),'_node_out.txt'];
FEB_struct.run_output_names={run_node_output_name};
FEB_struct.output_types={'node_data'};
FEB_struct.data_types={'ux;uy;uz'};

%Control section
FEB_struct.Control.AnalysisType='static';
FEB_struct.Control.Properties={'time_steps','step_size',...
    'max_refs','max_ups',...
    'dtol','etol','rtol','lstol'};
FEB_struct.Control.Values={10,0.1,...
    25,0,...
    0.001,0.01,0,0.9};
FEB_struct.Control.TimeStepperProperties={'dtmin','dtmax','max_retries','opt_iter','aggressiveness'};
FEB_struct.Control.TimeStepperValues={1e-5, 0.1, 5, 5, 1};

%Load curves
FEB_struct.LoadData.LoadCurves.id=1;
FEB_struct.LoadData.LoadCurves.type={'smooth'};
FEB_struct.LoadData.LoadCurves.loadPoints={[0 0;1 1]};

%% SAVING .FEB FILE

FEB_struct.disp_opt=0; %Display waitbars option
febStruct2febFile_v1p2(FEB_struct);

%% RUNNING FEBIO JOB

% FEBioRunStruct.FEBioPath='C:\Progra~1\FEBio1p8\febio.exe';
FEBioRunStruct.run_filename=FEB_struct.run_filename;
FEBioRunStruct.run_logname=FEB_struct.run_logname;
FEBioRunStruct.disp_on=1; 
FEBioRunStruct.disp_log_on=1; 
FEBioRunStruct.t_check=0.25; %Time for checking log file (dont set too small)
FEBioRunStruct.maxtpi=1e99; %Max analysis time
FEBioRunStruct.maxLogCheckTime=3;
%-------------------------------------------------------------------
[rundFlag]=runMonitorFEBio(FEBioRunStruct);%START FEBio NOW!!!!!!!!
%------------------------------------------------------------------

%% IMPORTING NODAL DISPLACEMENT RESULTS
% Importing nodal displacements from a log file
[~, N_disp_mat,~]=importFEBio_logfile(FEB_struct.run_output_names{1}); %Nodal displacements

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

hf1=figuremax(figColor,figColorDef);
title('Cut view of deformed model showing internal results','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;

hps=patch('Faces',Fs,'Vertices',VT_def,'FaceColor','flat','FaceVertexCData',Cs);

view(3); axis tight;  axis equal;  grid on;
colormap jet; colorbar; shading interp;
camlight headlight;
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

hf1=figuremax(figColor,figColorDef);
title('Outer surface only with distance metric','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;

hps=patch('Faces',F_free,'Vertices',VT_def,'FaceColor','flat','CData',CF);
hps=patch('Faces',F_free,'Vertices',VT,'FaceColor',0.5.*ones(1,3),'FaceAlpha',0.5,'EdgeColor','none');

view(3); axis tight;  axis equal;  grid on;
colormap jet; colorbar; caxis([0 max(CF(:))]);
camlight headlight;
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