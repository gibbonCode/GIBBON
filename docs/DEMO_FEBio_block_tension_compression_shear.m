%% DEMO_FEBio_block_tension_compression_shear
% Below is a demonstration for:
% 1) Building an FEBio model for tension, compression, and shear
% 2) Running the model
% 3) Importing and plotting displacement and force results

%%

clear; close all; clc;

%%
% Plot settings
fontSize=20;
faceAlpha1=0.5;
faceAlpha2=1;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize=40;
lineWidth=3;
plotColors=gjet(4);

%%
% Control parameters

% path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

modelNameEnd='tempModel';
modelName=fullfile(savePath,modelNameEnd);

%FEA control
numSteps=20;

%Specifying dimensions and number of elements
sampleWidth=36;
sampleThickness=sampleWidth; 
sampleHeight=sampleWidth;
pointSpacings=6*ones(1,3);
initialArea=sampleWidth*sampleThickness;

numElementsWidth=round(sampleWidth/pointSpacings(1));
numElementsThickness=round(sampleThickness/pointSpacings(2));
numElementsHeight=round(sampleHeight/pointSpacings(3));

stretchLoad=1.2;
displacementMagnitude=(stretchLoad*sampleHeight)-sampleHeight;

%Material parameter set
c1=1e-3; 
m1=6;
k=c1*20;

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
cFigure;
title('Model surfaces','FontSize',fontSize);
hold on;
gpatch(Fb,V,faceBoundaryMarker,'k',faceAlpha2,edgeWidth);

colormap(gjet(6)); icolorbar;
axisGeom(gca,fontSize);
camlight headlight;
drawnow;

%% DEFINE BC's

%Define supported node sets
logicFace=faceBoundaryMarker==5;
Fr=Fb(logicFace,:);
bcSupportList=unique(Fr(:));

%Prescribed displacement nodes
logicPrescribe=faceBoundaryMarker==6;
Fr=Fb(logicPrescribe,:);
bcPrescribeList=unique(Fr(:));

bcPrescribeMagnitudes=displacementMagnitude(ones(1,numel(bcPrescribeList)),:);

%%
% Visualize BC's
cFigure;
title('Boundary conditions','FontSize',fontSize);
hold on;

gpatch(Fb,V,0.5*ones(1,3),'k',faceAlpha1);
plotV(V(bcSupportList,:),'k.','MarkerSize',markerSize);
plotV(V(bcPrescribeList,:),'b.','MarkerSize',markerSize);

axisGeom(gca,fontSize);
camlight headlight;
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
% FEB_struct.Materials{1}.Type='Ogden';
% FEB_struct.Materials{1}.Name='Block_material';
% FEB_struct.Materials{1}.Properties={'c1','m1','c2','m2','k'};
% FEB_struct.Materials{1}.Values={c1,m1,c1,-m1,k};
FEB_struct.Materials{1}.Type='Ogden unconstrained';
FEB_struct.Materials{1}.Name='Block_material';
FEB_struct.Materials{1}.Properties={'c1','m1','c2','m2','cp'};
FEB_struct.Materials{1}.Values={c1,m1,c1,-m1,k};

%Step specific control sections
FEB_struct.Step{1}.Control.AnalysisType='static';
FEB_struct.Step{1}.Control.Properties={'time_steps','step_size',...
    'max_refs','max_ups',...
    'dtol','etol','rtol','lstol'};

FEB_struct.Step{1}.Control.Values={numSteps,1/numSteps,25,5,0.001,0.01,0,0.9};
FEB_struct.Step{1}.Control.TimeStepperProperties={'dtmin','dtmax','max_retries','opt_iter','aggressiveness'};
FEB_struct.Step{1}.Control.TimeStepperValues={(1/(100*numSteps)),1/numSteps,5,10,1};
FEB_struct.Step{2}=FEB_struct.Step{1};
FEB_struct.Step{3}=FEB_struct.Step{1};
FEB_struct.Step{4}=FEB_struct.Step{1};
FEB_struct.Step{5}=FEB_struct.Step{1};

%Defining node sets
FEB_struct.Geometry.NodeSet{1}.Set=bcSupportList;
FEB_struct.Geometry.NodeSet{1}.Name='bcSupportList';
FEB_struct.Geometry.NodeSet{2}.Set=bcPrescribeList;
FEB_struct.Geometry.NodeSet{2}.Name='bcPrescribeList';

%Adding BC information
FEB_struct.Boundary.Fix{1}.bc='x';
FEB_struct.Boundary.Fix{1}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
FEB_struct.Boundary.Fix{2}.bc='y';
FEB_struct.Boundary.Fix{2}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;
FEB_struct.Boundary.Fix{3}.bc='z';
FEB_struct.Boundary.Fix{3}.SetName=FEB_struct.Geometry.NodeSet{1}.Name;

%Prescribed BC's

%STEP 1 Tension
FEB_struct.Step{1}.Boundary.Prescribe{1}.Set=bcPrescribeList;
FEB_struct.Step{1}.Boundary.Prescribe{1}.bc='z';
FEB_struct.Step{1}.Boundary.Prescribe{1}.lc=1;
FEB_struct.Step{1}.Boundary.Prescribe{1}.nodeScale=displacementMagnitude(ones(numel(bcPrescribeList),1),1);
FEB_struct.Step{1}.Boundary.Prescribe{1}.Type='relative';

FEB_struct.Step{1}.Boundary.Prescribe{2}.Set=bcPrescribeList;
FEB_struct.Step{1}.Boundary.Prescribe{2}.bc='x';
FEB_struct.Step{1}.Boundary.Prescribe{2}.lc=1;
FEB_struct.Step{1}.Boundary.Prescribe{2}.nodeScale=zeros(numel(bcPrescribeList),1);
FEB_struct.Step{1}.Boundary.Prescribe{2}.Type='relative';

FEB_struct.Step{1}.Boundary.Prescribe{3}.Set=bcPrescribeList;
FEB_struct.Step{1}.Boundary.Prescribe{3}.bc='y';
FEB_struct.Step{1}.Boundary.Prescribe{3}.lc=1;
FEB_struct.Step{1}.Boundary.Prescribe{3}.nodeScale=zeros(numel(bcPrescribeList),1);
FEB_struct.Step{1}.Boundary.Prescribe{3}.Type='relative';

%STEP 2 Return form tension
FEB_struct.Step{2}.Boundary.Prescribe{1}.Set=bcPrescribeList;
FEB_struct.Step{2}.Boundary.Prescribe{1}.bc='z';
FEB_struct.Step{2}.Boundary.Prescribe{1}.lc=2;
FEB_struct.Step{2}.Boundary.Prescribe{1}.nodeScale=-displacementMagnitude(ones(numel(bcPrescribeList),1),1);
FEB_struct.Step{2}.Boundary.Prescribe{1}.Type='relative';

FEB_struct.Step{2}.Boundary.Prescribe{2}.Set=bcPrescribeList;
FEB_struct.Step{2}.Boundary.Prescribe{2}.bc='x';
FEB_struct.Step{2}.Boundary.Prescribe{2}.lc=2;
FEB_struct.Step{2}.Boundary.Prescribe{2}.nodeScale=zeros(numel(bcPrescribeList),1);
FEB_struct.Step{2}.Boundary.Prescribe{2}.Type='relative';

FEB_struct.Step{2}.Boundary.Prescribe{3}.Set=bcPrescribeList;
FEB_struct.Step{2}.Boundary.Prescribe{3}.bc='y';
FEB_struct.Step{2}.Boundary.Prescribe{3}.lc=2;
FEB_struct.Step{2}.Boundary.Prescribe{3}.nodeScale=zeros(numel(bcPrescribeList),1);
FEB_struct.Step{2}.Boundary.Prescribe{3}.Type='relative';

%STEP 3 Compression
FEB_struct.Step{3}.Boundary.Prescribe{1}.Set=bcPrescribeList;
FEB_struct.Step{3}.Boundary.Prescribe{1}.bc='z';
FEB_struct.Step{3}.Boundary.Prescribe{1}.lc=3;
FEB_struct.Step{3}.Boundary.Prescribe{1}.nodeScale=-displacementMagnitude(ones(numel(bcPrescribeList),1),1);
FEB_struct.Step{3}.Boundary.Prescribe{1}.Type='relative';

FEB_struct.Step{3}.Boundary.Prescribe{2}.Set=bcPrescribeList;
FEB_struct.Step{3}.Boundary.Prescribe{2}.bc='x';
FEB_struct.Step{3}.Boundary.Prescribe{2}.lc=3;
FEB_struct.Step{3}.Boundary.Prescribe{2}.nodeScale=zeros(numel(bcPrescribeList),1);
FEB_struct.Step{3}.Boundary.Prescribe{2}.Type='relative';

FEB_struct.Step{3}.Boundary.Prescribe{3}.Set=bcPrescribeList;
FEB_struct.Step{3}.Boundary.Prescribe{3}.bc='y';
FEB_struct.Step{3}.Boundary.Prescribe{3}.lc=3;
FEB_struct.Step{3}.Boundary.Prescribe{3}.nodeScale=zeros(numel(bcPrescribeList),1);
FEB_struct.Step{3}.Boundary.Prescribe{3}.Type='relative';

%STEP 4 Return from compression
FEB_struct.Step{4}.Boundary.Prescribe{1}.Set=bcPrescribeList;
FEB_struct.Step{4}.Boundary.Prescribe{1}.bc='z';
FEB_struct.Step{4}.Boundary.Prescribe{1}.lc=4;
FEB_struct.Step{4}.Boundary.Prescribe{1}.nodeScale=displacementMagnitude(ones(numel(bcPrescribeList),1),1);
FEB_struct.Step{4}.Boundary.Prescribe{1}.Type='relative';

FEB_struct.Step{4}.Boundary.Prescribe{2}.Set=bcPrescribeList;
FEB_struct.Step{4}.Boundary.Prescribe{2}.bc='x';
FEB_struct.Step{4}.Boundary.Prescribe{2}.lc=4;
FEB_struct.Step{4}.Boundary.Prescribe{2}.nodeScale=zeros(numel(bcPrescribeList),1);
FEB_struct.Step{4}.Boundary.Prescribe{2}.Type='relative';

FEB_struct.Step{4}.Boundary.Prescribe{3}.Set=bcPrescribeList;
FEB_struct.Step{4}.Boundary.Prescribe{3}.bc='y';
FEB_struct.Step{4}.Boundary.Prescribe{3}.lc=4;
FEB_struct.Step{4}.Boundary.Prescribe{3}.nodeScale=zeros(numel(bcPrescribeList),1);
FEB_struct.Step{4}.Boundary.Prescribe{3}.Type='relative';

%STEP 5 Shear
FEB_struct.Step{5}.Boundary.Prescribe{1}.Set=bcPrescribeList;
FEB_struct.Step{5}.Boundary.Prescribe{1}.bc='z';
FEB_struct.Step{5}.Boundary.Prescribe{1}.lc=5;
FEB_struct.Step{5}.Boundary.Prescribe{1}.nodeScale=zeros(numel(bcPrescribeList),1);
FEB_struct.Step{5}.Boundary.Prescribe{1}.Type='relative';

FEB_struct.Step{5}.Boundary.Prescribe{2}.Set=bcPrescribeList;
FEB_struct.Step{5}.Boundary.Prescribe{2}.bc='x';
FEB_struct.Step{5}.Boundary.Prescribe{2}.lc=5;
FEB_struct.Step{5}.Boundary.Prescribe{2}.nodeScale=displacementMagnitude(ones(numel(bcPrescribeList),1),1);
FEB_struct.Step{5}.Boundary.Prescribe{2}.Type='relative';

FEB_struct.Step{5}.Boundary.Prescribe{3}.Set=bcPrescribeList;
FEB_struct.Step{5}.Boundary.Prescribe{3}.bc='y';
FEB_struct.Step{5}.Boundary.Prescribe{3}.lc=5;
FEB_struct.Step{5}.Boundary.Prescribe{3}.nodeScale=zeros(numel(bcPrescribeList),1);
FEB_struct.Step{5}.Boundary.Prescribe{3}.Type='relative';


%Load curves
FEB_struct.LoadData.LoadCurves.id=[1 2 3 4 5];
FEB_struct.LoadData.LoadCurves.type={'linear','linear','linear','linear','linear'};
FEB_struct.LoadData.LoadCurves.loadPoints={[0 0;1 1];[0 0;1 0;2 1];[0 0;1 0;2 0;3 1];[0 0;1 0;2 0;3 0;4 1];[0 0;1 0;2 0;3 0;4 0; 5 1];};

%Adding output requests
FEB_struct.Output.VarTypes={'displacement','stress','relative volume'};

%Specify log file output
run_disp_output_name=[modelNameEnd,'_node_out.txt'];
run_force_output_name=[modelNameEnd,'_force_out.txt'];
FEB_struct.run_output_names={run_disp_output_name,run_force_output_name};
FEB_struct.output_types={'node_data','node_data'};
FEB_struct.data_types={'ux;uy;uz','Rx;Ry;Rz'};

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
    [time_mat, N_disp_mat,~]=importFEBio_logfile(fullfile(savePath,FEB_struct.run_output_names{1})); %Nodal displacements    
    
    %Remove nodal index column
    N_disp_mat=N_disp_mat(:,2:end,:);
    
    %Add initial state i.e. zero displacement
    sizImport=size(N_disp_mat); 
    sizImport(3)=sizImport(3)+1;
    N_disp_mat_n=zeros(sizImport);
    N_disp_mat_n(:,:,2:end)=N_disp_mat;
    N_disp_mat=N_disp_mat_n;
    
    %Add zero time point
    time_mat=[0; time_mat(:)]; %Time
    
    %% IMPORTING NODAL FORCES
    % Importing nodal forces from a log file
    [~, N_force_mat,~]=importFEBio_logfile(fullfile(savePath,FEB_struct.run_output_names{2})); %Nodal forces
        
    %Add initial state i.e. zero displacement
    sizImport=size(N_force_mat);
    sizImport(3)=sizImport(3)+1;
    N_force_mat_n=zeros(sizImport);
    N_force_mat_n(:,:,2:end)=N_force_mat;
    N_force_mat=N_force_mat_n;
    
    %% Plotting the deformed model
    
    DN_MAG=sqrt(sum(N_disp_mat.^2,2));
    DN=N_disp_mat(:,:,end);
    DN_magnitude=sqrt(sum(DN(:,3).^2,2));
    V_def=V+DN;
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
    X_DEF=V_DEF(:,1,:);
    Y_DEF=V_DEF(:,2,:);
    Z_DEF=V_DEF(:,3,:);
    [CF]=vertexToFaceMeasure(Fb,DN_magnitude);
    
    
    hf=cFigure;
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
    
    hp=gpatch(Fb,V_def,CF,'k',1);
    gpatch(Fb,V,0.5*ones(1,3),'k',0.25);    
    
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_MAG(:))]);
    axisGeom(gca,fontSize);
    axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);
    axis manual; 
    camlight headlight;   
    drawnow;
    
    animStruct.Time=time_mat;
    
    for qt=1:1:size(N_disp_mat,3)
        
        DN=N_disp_mat(:,:,qt);
        DN_magnitude=sqrt(sum(DN.^2,2));
        V_def=V+DN;
        [CF]=vertexToFaceMeasure(Fb,DN_magnitude);
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,CF}; %Property values for to set in order to animate
    end
        
    anim8(hf,animStruct);
    drawnow;
    
    %% Get force curves
    
    forceDataPrescribed=N_force_mat(bcPrescribeList,2:end,:);
   
    logicTension=time_mat>=0 & time_mat<=1;
    logicCompression=time_mat>=2 & time_mat<=3;
    logicShear=time_mat>=4 & time_mat<=5;
    
    forceTension=permute(sum(forceDataPrescribed(:,:,logicTension)),[3 2 1]);
    timeTension=time_mat(logicTension);
    timeTension=timeTension-timeTension(1);
    
    forceCompression=permute(sum(forceDataPrescribed(:,:,logicCompression)),[3 2 1]);
    timeCompression=time_mat(logicCompression);
    timeCompression=timeCompression-timeCompression(1);
    
    forceShear=permute(sum(forceDataPrescribed(:,:,logicShear)),[3 2 1]);
    timeShear=time_mat(logicShear);
    timeShear=timeShear-timeShear(1);
           
    %% Plotting force data
    
    cFigure; 
    ht=suptitle('Force history curves');
    ht.FontSize=fontSize;
    
    subplot(1,3,1); hold on;
    title('Tension');
    xlabel('Time [s]','FontSize',fontSize); ylabel('Force [N]','FontSize',fontSize); 
    h1=plot(timeTension,forceTension(:,1),'k--','Color',plotColors(1,:),'LineWidth',lineWidth);
    h2=plot(timeTension,forceTension(:,2),'k:','Color',plotColors(2,:),'LineWidth',lineWidth);
    h3=plot(timeTension,forceTension(:,3),'k-','Color',plotColors(4,:),'LineWidth',lineWidth);
    hl=legend([h1 h2 h3],'F_x','F_y','F_z','Location','southoutside','orientation','horizontal');
    set(gca,'FontSize',fontSize);
    axis square; axis tight; grid on;
    
    subplot(1,3,2); hold on;
    title('Compression')
    xlabel('Time [s]','FontSize',fontSize); ylabel('Force [N]','FontSize',fontSize); 
    h1=plot(timeCompression,forceCompression(:,1),'k--','Color',plotColors(1,:),'LineWidth',lineWidth);
    h2=plot(timeCompression,forceCompression(:,2),'k:','Color',plotColors(2,:),'LineWidth',lineWidth);
    h3=plot(timeCompression,forceCompression(:,3),'k-','Color',plotColors(4,:),'LineWidth',lineWidth);
    hl=legend([h1 h2 h3],'F_x','F_y','F_z','Location','southoutside','orientation','horizontal');
    set(gca,'FontSize',fontSize);
    axis square; axis tight; grid on;
    
    subplot(1,3,3); hold on;
    title('Shear')
    xlabel('Time [s]','FontSize',fontSize); ylabel('Force [N]','FontSize',fontSize); 
    h1=plot(timeShear,forceShear(:,1),'k--','Color',plotColors(1,:),'LineWidth',lineWidth);
    h2=plot(timeShear,forceShear(:,2),'k:','Color',plotColors(2,:),'LineWidth',lineWidth);
    h3=plot(timeShear,forceShear(:,3),'k-','Color',plotColors(4,:),'LineWidth',lineWidth);
    hl=legend([h1 h2 h3],'F_x','F_y','F_z','Location','southoutside','orientation','horizontal');
    set(gca,'FontSize',fontSize);
    axis square; axis tight; grid on;
    
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
% Copyright (C) 2017  Kevin Mattheus Moerman
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
