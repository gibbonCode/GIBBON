clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceAlpha1=0.3;
faceAlpha2=1;
cMap=gjet(4);
patchColor=cMap(1,:);
markerSize=25;

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement

%Path and file names for load curve data
loadNameLoadCurve=fullfile(defaultFolder,'data','fsi','Cornelissen_2018_data1.mat');

pointSpacing=4;
cylRadius=10;
cylHeight=cylRadius*7;
nt=ceil(cylHeight/pointSpacing);
numHeight=nt+iseven(nt);

distMarch=10; 

% Creating input structure
cylInputStruct.cylRadius=cylRadius;
cylInputStruct.numRadial=ceil((2*pi*cylRadius)/pointSpacing);
cylInputStruct.cylHeight=cylHeight;
cylInputStruct.numHeight=numHeight;
cylInputStruct.meshType='tri';
cylInputStruct.closeOpt=1;

layerThickness=cylRadius/5;
numElementsThickness=1;

pressure_offset=10000;

velocityData.min=0.1;
velocityData.max=0.5;
numWaves=1;

%%

%Material parameter set

%Bulk modulus factors
k_factor_vessel=1e2;
k_factor_clot=1e2;

%--> Vessel
c1_vessel=200000;%1000000; %Shear-modulus-like parameter
m1_vessel=2; %Material parameter setting degree of non-linearity
k_vessel=c1_vessel*k_factor_vessel; %Bulk modulus
ksi=c1_vessel/8;
alphaPar=2;
beta=2;

numMaterialsEnd=50;
enhancementFactor=100;
c1_vessel_range=linspace(c1_vessel,c1_vessel*enhancementFactor,numMaterialsEnd);
k_vessel_range=c1_vessel_range*k_factor_vessel;
materialIndexFluid=numMaterialsEnd+1;

%--> Fluid, fluid part
density_fluid  = 1060;
k_fluid        = 2.2e9;
mu0_Carreau    = 0.056;
mui_Carreau    = 0.00345;
lambda_Carreau = 3.313;
n_Carreau      = 0.3568;

%--> Fluid, solid part
density_neo=0;
E_neo=1e-9;
v_neo=0;

% FEA control settings
step_size_desired=0.004; %0.003;
numOutput=200;
max_refs=5;
max_ups=50;
diverge_reform=0;
reform_each_time_step=0;
dtol=0.001;
vtol=0.001;
ftol=0.001;
etol=0.01;
rtol=0.001;
lstol=0.9;
min_residual=1e-16;
max_residual=1e+10;
rhoi=0;
qnmethod=1;
max_retries=5;
opt_iter=53;
analysis_type='dynamic';
symmetric_stiffness=0;


%%
% Derive patch data for a cylinder
[F,V,C]=patchcylinder(cylInputStruct);
R=euler2DCM([0 0.5*pi 0]);
V=V*R;

V_regions=getInnerPoint(F,V); %Define region points
V_holes=[]; %Define hole points
[regionTetVolumes]=tetVolMeanEst(F,V); %Volume estimate for regular tets
stringOpt='-pq1.2AaY'; %Options for tetgen

%%
% Mesh using TetGen

%Create tetgen input structure
inputStruct.stringOpt=stringOpt; %Tetgen options
inputStruct.Faces=F; %Boundary faces
inputStruct.Nodes=V; %Nodes of boundary
inputStruct.faceBoundaryMarker=C;
inputStruct.regionPoints=V_regions; %Interior points for regions
inputStruct.holePoints=V_holes; %Interior points for holes
inputStruct.regionA=regionTetVolumes; %Desired tetrahedral volume for each region

% Mesh model using tetrahedral elements using tetGen
[meshOutput]=runTetGen(inputStruct); %Run tetGen

%%
% Access mesh output structure

E=meshOutput.elements; %The elements
V=meshOutput.nodes; %The vertices or nodes
CE=meshOutput.elementMaterialID; %Element material or region id
Fb=meshOutput.facesBoundary; %The boundary faces
Cb=meshOutput.boundaryMarker; %The boundary markers

%%
% Visualization

hf=cFigure;
subplot(1,2,1); hold on;
title('Input boundaries','FontSize',fontSize);
hp(1)=gpatch(Fb,V,Cb,'k',faceAlpha1);
hp(2)=plotV(V_regions,'r.','MarkerSize',markerSize);
legend(hp,{'Input mesh','Interior point(s)'},'Location','NorthWestOutside');
axisGeom(gca,fontSize); camlight headlight;
colormap(cMap); icolorbar;

hs=subplot(1,2,2); hold on;
title('Tetrahedral mesh','FontSize',fontSize);

% Visualizing using |meshView|
optionStruct.hFig=[hf,hs];
meshView(meshOutput,optionStruct);

axisGeom(gca,fontSize);
gdrawnow;

%%
% Visualizing meshed regions

cFigure;
gpatch(Fb,V,Cb);
colormap(gjet(3)); icolorbar;
axisGeom;
camlight headlight;
drawnow;

%% Create pentahedral skin
[Ep,Vp,Fq1,Fq2]=patchThick(Fb(Cb==1,:),V,-1,layerThickness,numElementsThickness);

%Use element2patch to get patch data
Fp=element2patch(Ep,[],'penta6');

%%
% Visualizing meshed regions

cFigure;
gpatch(Fb,V,'w','none',0.5);
gpatch(Fp,Vp,'rw','k',1);
axisGeom;
camlight headlight;
drawnow;

%%

F_inlet=Fb(Cb==2,:);
Eb_inlet=patchBoundary(F_inlet);
indBoundaryNodes_inlet=edgeListToCurve(Eb_inlet);
F_outlet=Fb(Cb==3,:);
Eb_outlet=patchBoundary(F_outlet);
indBoundaryNodes_outlet=edgeListToCurve(Eb_outlet);

%%
% Visualizing meshed regions

cFigure; hold on;
gpatch(Fb,V,'w','none',0.5);
hp1=gpatch(F_inlet,V,'rw','k',1);
hp2=gpatch(F_outlet,V,'bw','k',1);
hp3=plotV(V(indBoundaryNodes_inlet,:),'r-','LineWidth',3);
hp4=plotV(V(indBoundaryNodes_outlet,:),'b-','LineWidth',3);
legend([hp1 hp2 hp3 hp4],{'Inlet faces','Outlet faces','Inlet boundary nodes','Outlet boundary nodes'})
axisGeom;
camlight headlight;
drawnow;

%% set-up load curves

% Load load-curve data
loadCurveData=load(loadNameLoadCurve);
t=loadCurveData.x; %Time
a=loadCurveData.y; %Amplitude

tDevelopStart=(max(distMarch)./velocityData.min);
tWaitStart=0.5*tDevelopStart;
tWaitIntermediate=(max(t(:))-min(t(:)))/6;

sigmoidRampDurationVelocity=tDevelopStart/2;
sigmoidRampDuration_P0=sigmoidRampDurationVelocity+(tDevelopStart-sigmoidRampDurationVelocity)/2;

w=velocityData.max-velocityData.min;

%in-vivo part

a=a-min(a(:));
a=a./max(a(:));
a=a.*w;
a=a+velocityData.min;

%Sigmoid part
s=10;
T=linspace(-sigmoidRampDurationVelocity/2,sigmoidRampDurationVelocity/2,100);

S=tanh(T.*s);
S=S-min(S(:));
S=S./max(S(:));
S=S.*velocityData.min;
S(end+1)=S(end);
T=T-T(1);
T(end+1)=T(end)+tWaitStart;

t=t+T(end);
t(end+1)=t(end)+tWaitIntermediate;
a(end+1)=a(end);

TT=[T(:);];
AA=[S(:);];

for q=1:1:numWaves
    TT=[TT;t(:)-min(t(:))+max(TT(:))];
    AA=[AA;a(:)];
end

[TT,indUni]=unique(TT);
AA=AA(indUni);
loadCurve_velocity=[TT(:) AA(:)];

dt=mean(diff(loadCurve_velocity(:,1)));
n=ceil(max(loadCurve_velocity(:,1))/dt);
loadCurve_velocity=evenlySampleCurve(loadCurve_velocity,n,'pchip',0);

%Sigmoid part
s=15;
T=linspace(-sigmoidRampDuration_P0/2,sigmoidRampDuration_P0/2,100);

S=tanh(T.*s);
S=S-min(S(:));
S=S./max(S(:));
%     S=S.*velocityData.min;
%     S(end+1)=S(end);
T=T-T(1);
%     T(end+1)=T(end)+tWaitStart;

dt=mean(diff(loadCurve_velocity(:,1)));
n=ceil(max(T(:))/dt);
loadCurve_other=[T(:) S(:)];
loadCurve_other=evenlySampleCurve(loadCurve_other,n,'pchip',0);

%%

cFigure;
subplot(1,3,1);hold on;
plot(loadCurve_velocity(:,1),loadCurve_velocity(:,2),'r.-','MarkerSize',15,'LineWidth',2);
set(gca,'FontSize',35);
grid on; box on; axis square;
set(gca,'XTick',0:0.5:max(loadCurve_velocity(:,1)));

subplot(1,3,2);hold on;
plot(loadCurve_other(:,1),loadCurve_other(:,2),'g.-','MarkerSize',15,'LineWidth',2);
set(gca,'FontSize',35);
grid on; box on; axis square;
set(gca,'XTick',0:0.1:max(loadCurve_other(:,1)));

subplot(1,3,3);hold on;
plot(loadCurve_velocity(:,1),loadCurve_velocity(:,2),'r.-','MarkerSize',15,'LineWidth',2);
plot(loadCurve_other(:,1),loadCurve_other(:,2).*velocityData.min,'g.-','MarkerSize',15,'LineWidth',2);
set(gca,'FontSize',35);
grid on; box on; axis square;
set(gca,'XTick',0:0.5:max(loadCurve_velocity(:,1)));

drawnow;

%%

cFigure; hold on;
title('Inlet fluid velocity profile');
plot(loadCurve_velocity(:,1),loadCurve_velocity(:,2),'r-','LineWidth',5);
set(gca,'FontSize',35);
grid on; box on; %axis square;
axis tight;
set(gca,'XTick',0:0.25:max(loadCurve_velocity(:,1)));
set(gca,'YTick',round(100*linspace(0,max(loadCurve_velocity(:,2)),11))/100);
ylim([0 0.5]);
xlabel('Time [s]');
ylabel('Velocity [m/s]');
ha=gca;
ha.LineWidth=2;
drawnow;

dasfasas

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings
[febio_spec]=febioStructTemplate;

%febio_spec version
febio_spec.ATTR.version='3.0';

%Module section
febio_spec.Module.ATTR.type='solid';

%Control section
febio_spec.Control.analysis='STATIC';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.solver.max_ups=max_ups;
febio_spec.Control.solver.symmetric_stiffness=symmetric_stiffness;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax;
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;

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
febio_spec.Material.material{2}.ATTR.type='rigid body';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.density=1;
febio_spec.Material.material{2}.center_of_mass=center_of_mass;

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
febio_spec.Mesh.Elements{2}.ATTR.type='tri3'; %Element type
febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=E2; %The element matrix

% -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.node.ATTR.id=bcSupportList(:);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat=materialName1;

febio_spec.MeshDomains.SolidDomain{2}.ATTR.name=partName2;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat=materialName2;

%Boundary condition section
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.type='fix';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.dofs='x,y,z';


%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0; 1 1];

febio_spec.LoadData.load_controller{2}.ATTR.id=2;
febio_spec.LoadData.load_controller{2}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{2}.points.point.VAL=[0 0; 1 1];

febio_spec.LoadData.load_controller{3}.ATTR.id=3;
febio_spec.LoadData.load_controller{3}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{3}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{3}.points.point.VAL=[0 0; 1 1];


%Output section
% -> log file

% febio_spec.Output.logfile.ATTR.file=febioLogFileName;
% febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
% febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
% febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

febio_spec.Output.plotfile.var{end+1}.ATTR.type='velocity';
febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid pressure';
febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid density';
febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid dilatation';
febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid stress';
febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid velocity';
febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid volume ratio';
febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid shear viscosity';
febio_spec.Output.plotfile.var{end+1}.ATTR.type='fluid vorticity';
febio_spec.Output.plotfile.var{end+1}.ATTR.type='nodal fluid velocity';
febio_spec.Output.plotfile.var{end+1}.ATTR.type='nodal relative fluid velocity';
febio_spec.Output.plotfile.var{end+1}.ATTR.type='relative fluid velocity';

%% Quick viewing of the FEBio input file structure
% The |febView| function can be used to view the xml structure in a MATLAB
% figure window.

%%
% |febView(febio_spec); %Viewing the febio file|

%% Exporting the FEBio input file
% Exporting the febio_spec structure to an FEBio input file is done using
% the |febioStruct2xml| function.

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode


% %% Running the FEBio analysis
% % To run the analysis defined by the created FEBio input file the
% % |runMonitorFEBio| function is used. The input for this function is a
% % structure defining job settings e.g. the FEBio input file name. The
% % optional output runFlag informs the user if the analysis was run
% % succesfully.
%
% febioAnalysis.run_filename=febioFebFileName; %The input file name
% febioAnalysis.run_logname=febioLogFileName; %The name for the log file
% febioAnalysis.disp_on=1; %Display information on the command window
% febioAnalysis.runMode='internal';%'internal';
%
% [runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

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
