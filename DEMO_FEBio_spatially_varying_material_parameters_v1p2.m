%% DEMO_spatially_varying_material_parameters
% This demonstration shows how to generate a model with spatially varying
% material parameters.

%%
clear; close all; clc; 

%%
% Plot settings
figColor='w'; figColorDef='white';
fontSize=15;
faceAlpha1=1;
faceAlpha2=0.5;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize1=25;
boneColor=[1 0.9 0.8];
softColor=[0.6 0.08 0.08];

cMap=linspacen(softColor(:),boneColor(:),250)';

%%
% path names
filePath=mfilename('fullpath');
savePath=fullfile(fileparts(filePath),'data','temp');
modelName=fullfile(savePath,'tempModel');

%% DEFINING AND VISUALIZING THE PARAMETER MAP
% A trabecular structure is here simulated using the "Gyroid triply
% periodic surface" function.  

%Define Mooney-Rivlin parameter
nPar=15; %Number of parameter levels
minC=1e-5; %minimum value
maxC=1e-3; %Maximum value
c1_range=linspace(minC,maxC,nPar); %Value range

n=20;
[X,Y,Z]=meshgrid(linspace(-2*pi,2*pi,n));
V=[X(:) Y(:) Z(:)];
[R,~]=euler2DCM([0.25*pi 0.25*pi 0.25*pi]);
V=(R*V')';
S=triplyPeriodicMinimal(V(:,1),V(:,2),V(:,3),'g');
S=reshape(S,size(X));

%% 
% VISUALIZING THE MAPPING

[F,V,C]=ind2patch(true(size(S)),S,'vb'); 
[C_rgb]=gray2RGBColorMap(C,jet(250),[min(S(:)) max(S(:))]);

[Fs1,Vs1,Cs1]=ind2patch(S>=0.6,S,'vb'); 
[Fs2,Vs2,Cs2]=ind2patch(S<0.6,S,'vb'); 

figuremax(figColor,figColorDef);

subplot(1,2,1);
title('Stiff network','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
patch('Faces',Fs1,'Vertices',Vs1,'FaceColor','flat','CData',Cs1,'EdgeColor','k','FaceAlpha',1);
axis equal; view(3); axis tight; axis vis3d; grid on; view([-20,22]);
colormap(cMap); caxis([min(S(:)) max(S(:))]); colorbar; 
set(gca,'FontSize',fontSize);

subplot(1,2,2);
title('Soft matrix','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
patch('Faces',Fs2,'Vertices',Vs2,'FaceColor','flat','CData',Cs2,'EdgeColor','k','FaceAlpha',1);
axis equal; view(3); axis tight; axis vis3d; grid on; view([-20,22]);
colormap(cMap); caxis([min(S(:)) max(S(:))]); colorbar;
set(gca,'FontSize',fontSize);

drawnow;

%% BUILD MODEL

%Create hexahedral elements with function based colors
[E,V,C]=ind2patch(true(size(S)),S,'hu'); 

%Define element parameter mapping
elementMaterialIndices=C;
elementMaterialIndices=elementMaterialIndices-min(elementMaterialIndices(:));
elementMaterialIndices=elementMaterialIndices./max(elementMaterialIndices(:));
elementMaterialIndices=round(elementMaterialIndices.*(nPar-1))+1;

[F,PF]=element2patch(E,elementMaterialIndices);

logicTopNodes=abs(V(:,3)-max(V(:,3)))<=max(eps(V(:,3)));

logicBottomNodes=abs(V(:,3)-min(V(:,3)))<=max(eps(V(:,3)));


figuremax(figColor,figColorDef);
title('Gyroid derived model of trabecular structure','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',PF,'edgeColor',edgeColor,'lineWidth',edgeWidth,'FaceAlpha',1);
plotV(V(logicTopNodes,:),'k.','MarkerSize',markerSize1);
plotV(V(logicBottomNodes,:),'k.','MarkerSize',markerSize1);
colormap(cMap); caxis([min(elementMaterialIndices(:)) max(elementMaterialIndices(:))]); colorbar;

axis equal; view(3); axis tight; axis vis3d; grid on;  
set(gca,'FontSize',fontSize);

%% SET UP BOUNDARY CONDITIONS

%List of nodes to fix
bcFixList=find(logicBottomNodes);

%List of nodes to prescribe displacement for
bcPrescribeList=find(logicTopNodes);

%Define displacement magnitudes
bcPrescribedMagnitudes=zeros(numel(bcPrescribeList),1);
bcPrescribedMagnitudes(:,3)=3;

%% CONSTRUCTING FEB MODEL

% Defining file names
FEB_struct.run_filename=[modelName,'.feb']; %FEB file name
FEB_struct.run_logname=[modelName,'.txt']; %FEBio log file name

%Creating FEB_struct
FEB_struct.Geometry.Nodes=V;
FEB_struct.Geometry.Elements={E}; %The element sets
FEB_struct.Geometry.ElementType={'hex8'}; %The element types
FEB_struct.Geometry.ElementMat={elementMaterialIndices};

% DEFINING SPATIALLY VARYING MATERIAL SET

%Entries not varying per element
%Material 1
k_factor=100;
MatQ.type='Mooney-Rivlin';
MatQ.props={'c1','c2','k'};
MatQ.aniso_type='none';

for q=1:1:nPar      
    %Defining material parameters
    c1=c1_range(q);
    k=c1*k_factor;
    MatQ.vals={c1,0,k};
    FEB_struct.Materials{q}=MatQ;
end

%Adding BC information
FEB_struct.Boundary.FixList={bcFixList};
FEB_struct.Boundary.FixType={'xyz'};

FEB_struct.Boundary.PrescribeList={bcPrescribeList,bcPrescribeList,bcPrescribeList};
FEB_struct.Boundary.PrescribeType={'x','y','z'};

FEB_struct.Boundary.PrescribeValues={bcPrescribedMagnitudes(:,1),bcPrescribedMagnitudes(:,2),bcPrescribedMagnitudes(:,3)};
FEB_struct.Boundary.LoadCurveIds=[1 1 1];

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

FEB_struct.disp_opt=0; %Display waitbars option

%% SAVING .FEB FILE

febStruct2febFile_v1p2(FEB_struct);

%% RUNNING FEBIO JOB

% FEBioRunStruct.FEBioPath='C:\Progra~1\FEBio1p8\febio.exe';
FEBioRunStruct.run_filename=FEB_struct.run_filename;
FEBioRunStruct.run_logname=FEB_struct.run_logname;
FEBioRunStruct.disp_on=1; 
FEBioRunStruct.disp_log_on=1; 
FEBioRunStruct.t_check=0.25; %Time for checking log file (dont set too small)
FEBioRunStruct.maxtpi=1e99; %Max analysis time
FEBioRunStruct.maxLogCheckTime=3; %Max log file checking time

%-------------------------------------------------------------------
[rundFlag]=runMonitorFEBio(FEBioRunStruct);%START FEBio NOW!!!!!!!!
%------------------------------------------------------------------

%% IMPORTING NODAL DISPLACEMENT RESULTS
% Importing nodal displacements from a log file
[~, N_disp_mat,~]=importFEBio_logfile(FEB_struct.run_output_names{1}); %Nodal displacements

DN=N_disp_mat(:,2:end,end); %Final nodal displacements

%% CREATING NODE SET IN DEFORMED STATE
V_def=V+DN;
DN_magnitude=sqrt(sum(DN.^2,2));

%%
% Plotting the deformed model

[CF]=vertexToFaceMeasure(F,DN_magnitude);

hf1=figuremax(figColor,figColorDef);
title('The deformed model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;

hps=patch('Faces',F,'Vertices',V_def,'FaceColor','flat','CData',CF,'lineWidth',edgeWidth,'edgeColor',edgeColor,'FaceAlpha',faceAlpha1);

view(3); axis tight;  axis equal;  grid on;
colormap jet; colorbar;
% camlight headlight;
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