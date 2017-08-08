%% DEMO_FEBio_hexLattice_compression
% Below is a demonstration for: 
%
% * The creation of an FEBio model whereby force is applied to a selection
% of nodes, in this case to the end of a bar
% *  Running an FEBio job with MATLAB
% *  Importing FEBio results into MATLAB

%%

clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceAlpha1=0.8;
faceAlpha2=1;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize=50;
cMap=gjet(4);

%%
% Control parameters

% path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

modelNameEnd='tempModel';
modelName=fullfile(savePath,modelNameEnd);

%Specifying dimensions and number of elements
sampleSize=10;
stretchLoad=0.8;
displacementMagnitude=[0 0 (stretchLoad*sampleSize)-sampleSize];

%Material parameters
c=1e-3;
k=1e3*c;
m=2;

%Analysis contol parameters
nSteps=15;
timeStep=1/nSteps;
max_refs=25;
max_ups=0;

%%
% Creating example geometry. 
X=[-1;  1; 1; -1; -1;  1; 1; -1;];
Y=[-1; -1; 1;  1; -1; -1; 1;  1;];
Z=[-1; -1;-1; -1;  1;  1; 1;  1;];
V=sampleSize*[X(:) Y(:) Z(:)]/2;
E=1:8; %Element description of the 8-node cube (hexahedral element)

[E,V,C]=subHex(E,V); %Subdevide into 8 sub-cubes
[E,V,C]=hex2tet(E,V,C,1); %Convert to tetrahedral elements
[F,~]=element2patch(E,C); %Patch data for plotting

%%
% Create lattice structure
controlParameter.growSteps=0; %0 is normal, positive or negative integers increase or decrease the edge lattice thickness respectively
controlParameter.latticeSide=2; %Empty outputs both, 1=side 1 the edge lattice, 2=side 2 the dual lattice to the edge lattice
[Es,Vs,Cs]=element2HexLattice(E,V,controlParameter); %Get lattice structure
elementMaterialIndices=ones(size(Es,1),1);

% Create patch Data for visualization
[Fs,CsF]=element2patch(Es,Cs); %Patch data for plotting

indB=tesBoundary(Fs,Vs);
Fb=Fs(indB,:);

%%
% Visualizing input mesh and lattic structures

cFigure;
hs=subplot(1,2,1); 
title('The input mesh','fontSize',fontSize)
hold on;
gpatch(F,V,0.5*ones(1,3),'k',0.5);
axisGeom(gca,fontSize); 
camlight headlight; lighting flat;

subplot(1,2,2); 
title('Lattice side 1','fontSize',fontSize)
hold on;
% gpatch(Fs,Vs,'b');
gpatch(Fb,Vs,'b');
% patchNormPlot(Fs,Vs);
axisGeom(gca,fontSize); 
camlight headlight; lighting flat;

drawnow;

%%

faceBoundaryMarker=zeros(size(Fb,1),1);

%Find top and bottom face sets
[Nb]=patchNormal(Fb,Vs);
zVec=[0 0 1];
d=dot(Nb,zVec(ones(size(Nb,1),1),:),2);
Z=Vs(:,3);
ZF=mean(Z(Fb),2);
logicTop_Fb=(d>0.9) & ZF>=(max(Vs(:,3))-eps(sampleSize));
logicBottom_Fb=(d<-0.9) & ZF<=(min(Vs(:,3))+eps(sampleSize));

xVec=[1 0 0];
d=dot(Nb,xVec(ones(size(Nb,1),1),:),2);
X=Vs(:,1);
XF=mean(X(Fb),2);
logicSides_Fb1=(d>0.9) & XF>=(max(Vs(:,1))-eps(sampleSize));
logicSides_Fb2=(d<-0.9) & XF<=(min(Vs(:,1))+eps(sampleSize));

yVec=[0 1 0];
d=dot(Nb,yVec(ones(size(Nb,1),1),:),2);
Y=Vs(:,2);
YF=mean(Y(Fb),2);
logicSides_Fb3=(d>0.9) & YF>=(max(Vs(:,2))-eps(sampleSize));
logicSides_Fb4=(d<-0.9) & YF<=(min(Vs(:,2))+eps(sampleSize));

logicSides_Fb=logicSides_Fb1 | logicSides_Fb2 | logicSides_Fb3 | logicSides_Fb4;

faceBoundaryMarker(logicBottom_Fb)=1;
faceBoundaryMarker(logicTop_Fb)=2;
faceBoundaryMarker(logicSides_Fb)=3;

cFigure;
hold on;
gpatch(Fs,Vs,0.5*ones(1,3),'none',0.2);
gpatch(Fb(faceBoundaryMarker==1,:),Vs,'r','k',1);
gpatch(Fb(faceBoundaryMarker==2,:),Vs,'b','k',1);
gpatch(Fb(faceBoundaryMarker==3,:),Vs,'g','k',1);
axisGeom(gca,fontSize); 
colormap(cMap);
camlight headlight; lighting flat;
drawnow;

%% DEFINE BC's

%Supported nodes
logicRigid=faceBoundaryMarker==1;
Fr=Fb(logicRigid,:);
bcRigidList=unique(Fr(:));

%Prescribed force nodes
logicPrescribe=faceBoundaryMarker==2;
Fr=Fb(logicPrescribe,:);
bcPrescribeList=unique(Fr(:));
bcPrescribeMagnitudes=displacementMagnitude(ones(1,numel(bcPrescribeList)),:);

%%
% Visualize BC's

cFigure; hold on;
title('Boundary conditions','FontSize',fontSize);
gpatch(Fs,Vs,0.5*ones(1,3),'k',0.4);
plotV(Vs(bcRigidList,:),'b.','MarkerSize',markerSize);
plotV(Vs(bcPrescribeList,:),'r.','MarkerSize',markerSize);
axisGeom;
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow; 

%%

[Fb_clean,Vb_clean,indFix]=patchCleanUnused(Fb,Vs);

cPar.Method='HC';
cPar.n=15;
indKeep=Fb_clean(faceBoundaryMarker>0,:);
indKeep=unique(indKeep(:));
cPar.RigidConstraints=indKeep;

[Vb_clean]=tesSmooth(Fb_clean,Vb_clean,[],cPar);
ind=Fb(:);
ind=unique(ind(:));
Vs(ind,:)=Vb_clean;

cFigure; hold on;
title('Boundary conditions','FontSize',fontSize);
gpatch(Fs,Vs,'b','k',1);
axisGeom;
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow; 

%% CONSTRUCTING FEB MODEL

FEB_struct.febio_spec.version='2.0';
FEB_struct.Module.Type='solid';

% Defining file names
FEB_struct.run_filename=[modelName,'.feb']; %FEB file name
FEB_struct.run_logname=[modelName,'.txt']; %FEBio log file name

%Geometry section
FEB_struct.Geometry.Nodes=Vs;
FEB_struct.Geometry.Elements={Es}; %The element sets
FEB_struct.Geometry.ElementType={'hex8'}; %The element types
FEB_struct.Geometry.ElementMat={elementMaterialIndices};
FEB_struct.Geometry.ElementsPartName={'Bar'};

%Material section
FEB_struct.Materials{1}.Type='Ogden';
FEB_struct.Materials{1}.Properties={'c1','c2','m1','m2','k'};
FEB_struct.Materials{1}.Values={c,c,m,-m,k};

%Step specific control sections
FEB_struct.Control.AnalysisType='static';
FEB_struct.Control.Properties={'time_steps','step_size',...
    'max_refs','max_ups',...
    'dtol','etol','rtol','lstol'};

FEB_struct.Control.Values={nSteps,timeStep,...
    max_refs,max_ups,...
    0.001,0.01,0,0.9};
FEB_struct.Control.TimeStepperProperties={'dtmin','dtmax','max_retries','opt_iter'};
FEB_struct.Control.TimeStepperValues={timeStep/100,timeStep, 5, 5};

%Defining node sets
FEB_struct.Geometry.NodeSet{1}.Set=bcRigidList;
FEB_struct.Geometry.NodeSet{1}.Name='bcRigidList';
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
FEB_struct.Boundary.Prescribe{1}.Set=bcPrescribeList;
FEB_struct.Boundary.Prescribe{1}.bc='x';
FEB_struct.Boundary.Prescribe{1}.lc=1;
FEB_struct.Boundary.Prescribe{1}.nodeScale=bcPrescribeMagnitudes(:,1);
FEB_struct.Boundary.Prescribe{1}.Type='relative';

FEB_struct.Boundary.Prescribe{2}.Set=bcPrescribeList;
FEB_struct.Boundary.Prescribe{2}.bc='y';
FEB_struct.Boundary.Prescribe{2}.lc=1;
FEB_struct.Boundary.Prescribe{2}.nodeScale=bcPrescribeMagnitudes(:,2);
FEB_struct.Boundary.Prescribe{2}.Type='relative';

FEB_struct.Boundary.Prescribe{3}.Set=bcPrescribeList;
FEB_struct.Boundary.Prescribe{3}.bc='z';
FEB_struct.Boundary.Prescribe{3}.lc=1;
FEB_struct.Boundary.Prescribe{3}.nodeScale=bcPrescribeMagnitudes(:,3);
FEB_struct.Boundary.Prescribe{3}.Type='relative';

%Load curves
FEB_struct.LoadData.LoadCurves.id=1;
FEB_struct.LoadData.LoadCurves.type={'linear'};
FEB_struct.LoadData.LoadCurves.loadPoints={[0 0;1 1;]};

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
    [~, N_disp_mat,~]=importFEBio_logfile(fullfile(savePath,FEB_struct.run_output_names{1})); %Nodal displacements    
    
    %% IMPORTING NODAL FORCES
    % Importing nodal forces from a log file
    [time_mat, N_force_mat,~]=importFEBio_logfile(fullfile(savePath,FEB_struct.run_output_names{2})); %Nodal forces
    time_mat=[0; time_mat(:)]; %Time
    
    N_disp_mat=N_disp_mat(:,2:end,:);
    sizImport=size(N_disp_mat);
    sizImport(3)=sizImport(3)+1;
    N_disp_mat_n=zeros(sizImport);
    N_disp_mat_n(:,:,2:end)=N_disp_mat;
    N_disp_mat=N_disp_mat_n;
    DN=N_disp_mat(:,:,end);
    DN_magnitude=sqrt(sum(DN(:,3).^2,2));
    Vs_def=Vs+DN;
    [CF]=vertexToFaceMeasure(Fb,DN_magnitude);
    
    %% Axis limits for plotting
    
    minD=min(Vs+min(N_disp_mat,[],3),[],1);
    maxD=max(Vs+max(N_disp_mat,[],3),[],1);
    axisLim=[minD(1) maxD(1) minD(2) maxD(2) minD(3) maxD(3)];
    
    %% Plotting the deformed model
    
    hf=cFigure;
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
    
    hp=gpatch(Fb,Vs_def,CF,'k',1);
%     gpatch(Fb,Vs,0.5*ones(1,3),'none',0.2);
    
    view(3); axis tight;  axis equal;  grid on; box on;
    colormap(gjet(250)); colorbar;
%     caxis([0 max(DN_magnitude)]);
%     view(130,25);
    camlight headlight; lighting flat;
    set(gca,'FontSize',fontSize);
    axis(axisLim);
    axis off;
    drawnow;

    animStruct.Time=time_mat;
    
    for qt=1:1:size(N_disp_mat,3)
        
        DN=N_disp_mat(:,:,qt);
        DN_magnitude=sqrt(sum(DN(:,3).^2,2));
        Vs_def=Vs+DN;
        [CF]=vertexToFaceMeasure(Fb,DN_magnitude);
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={Vs_def,CF}; %Property values for to set in order to animate
    end
        
    anim8(hf,animStruct);
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
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
