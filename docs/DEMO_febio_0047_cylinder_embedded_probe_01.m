%% DEMO_febio_0047_cylinder_embedded_probe_01
% Below is a demonstration for:
% 
% * Building geometry for a tissue segment with an embedded probe
% * Defining the boundary conditions 
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement results

%% Keywords
%
% * febio_spec version 2.5
% * febio, FEBio
% * probe
% * contact, sliding, sticky, friction
% * rigid body constraints
% * hexahedral elements, hex8
% * triangular elements, tri3
% * slab, block, rectangular
% * sphere
% * static, solid
% * hyperelastic, Ogden
% * displacement logfile
% * stress logfile

%%
clear; close all; clc;

%% 
% Plot settings
fontSize=15;
faceAlpha=1;
lineWidth1=1.5;
lineWidth2=3; 
markerSize1=15; 
markerSize2=30; 
edgeWidth=2; 
edgeColor='k';
faceAlpha1=1; 

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
febioLogFileName_strainEnergy=[febioFebFileNamePart,'_energy_out.txt']; %Log file name for exporting strain energy density
febioLogFileName_strain=[febioFebFileNamePart,'_strain_out.txt']; %Log file name for exporting strain

probeHeight=75;

probeRadius=2; % The radius of the hemi-spher portion
nRefine=1;  % Number of |subtri| refinements for icosahedron

pointSpacing=2; 
dAdd=7;
tissueRadius=probeRadius+dAdd;
tissueHeight=probeHeight+dAdd;
volumeFactor=5;

displacementMagnitude=-1;

%Material parameter set
c1=1e-3; %Shear-modulus-like parameter
m1=8; %Material parameter setting degree of non-linearity
k_factor=1e2; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus

% FEA control settings
numTimeSteps=20; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

%% Build probe

probeMeshInputStruct.sphereRadius=probeRadius;% => The radius of the hemi-spher portion
probeMeshInputStruct.nRefine=nRefine;% => Number of |subtri| refinements for icosahedron
probeMeshInputStruct.cylinderHeight=probeHeight-probeRadius;% => height of the cylinder part
probeMeshInputStruct.cylinderStepSize=[];% => Aproximate node spacing for cylinder portion
probeMeshInputStruct.patchType='tri';

[Fp,Vp,Cp]=hemiSphereCylMesh(probeMeshInputStruct);
Fp=fliplr(Fp); %Invert face orientation

Vp(:,3)=Vp(:,3)-max(Vp(:,3));

%Get top curve
Eb=patchBoundary(Fp,Vp);
indProbeTop=edgeListToCurve(Eb);
indProbeTop=indProbeTop(1:end-1);
Vst=Vp(indProbeTop,:);

%%
cFigure; hold on;
gpatch(Fp,Vp,'gw','k');
plotV(Vst,'b.-','lineWidth',lineWidth1,'MarkerSize',markerSize1);
patchNormPlot(Fp,Vp);
axisGeom(gca,fontSize);
drawnow; 

%% Build gel

%Sketching profile
ns=150;
t=linspace(0,2*pi,ns);
t=t(1:end-1);

x=tissueRadius*cos(t);
y=tissueRadius*sin(t);
z=zeros(size(x));
Vc=[x(:) y(:) z(:)];
np=ceil(max(pathLength(Vc))./pointSpacing);
[Vc]=evenlySampleCurve(Vc,np,'pchip',1);
 
% Extruding model
cPar.numSteps=round(tissueHeight/pointSpacing);
cPar.depth=tissueHeight; 
cPar.patchType='tri'; 
cPar.dir=-1;
cPar.closeLoopOpt=1; 
[Fg,Vg]=polyExtrude(Vc,cPar);
Fg=fliplr(Fg);

Vgb=Vg(cPar.numSteps:cPar.numSteps:end,:); 

Vgt=Vg(1:cPar.numSteps:end,:); 

%% Cap ends

regionCell={Vgt(:,[1 2]),Vst(:,[1 2])};
[Ft,Vt]=regionTriMesh2D(regionCell,pointSpacing,0,0);
Vt(:,3)=mean(Vgt(:,3));

regionCell={Vgb(:,[1 2])};
[Fb,Vb]=regionTriMesh2D(regionCell,pointSpacing,0,0);
Fb=fliplr(Fb); %flip face orientation
Vb(:,3)=mean(Vgb(:,3));

%%
% Visualize

cFigure; hold on;

gpatch(Fp,Vp,'rw','k',0.5);
gpatch(Fg,Vg,'gw','k',0.5);
gpatch(Fb,Vb,'bw','k',0.5);
gpatch(Ft,Vt,'bw','k',0.5);

plotV(Vgb,'b.-','lineWidth',lineWidth1,'MarkerSize',markerSize1);
plotV(Vgt,'b.-','lineWidth',lineWidth1,'MarkerSize',markerSize1);
plotV(Vst,'b.-','lineWidth',lineWidth1,'MarkerSize',markerSize1);

axisGeom(gca,fontSize);
drawnow; 

%% Merge model components

[F,V,C]=joinElementSets({Fg,Ft,Fb,Fp},{Vg,Vt,Vb,Vp});
[F,V]=mergeVertices(F,V);

%%

cFigure; 
subplot(1,2,1); hold on;
gpatch(F,V,C,'none',0.5);
axisGeom(gca,fontSize);
colormap gjet; icolorbar;

subplot(1,2,2); hold on;
gpatch(F,V,C);
patchNormPlot(F,V,2);
plotV(Vst,'b.-','lineWidth',lineWidth1,'MarkerSize',markerSize1);
axisGeom(gca,fontSize);
colormap gjet; icolorbar;

drawnow; 

%% Mesh solid using tetgen

%%
% Create tetgen meshing input structure

[regionA]=tetVolMeanEst(F,V); %Volume for a regular tet based on edge lengths
V_inner=getInnerPoint(F,V); %Interior point for region

inputStruct.stringOpt='-pq1.2AaY';
inputStruct.Faces=F;
inputStruct.Nodes=V;
inputStruct.holePoints=[];
inputStruct.faceBoundaryMarker=C; %Face boundary markers
inputStruct.regionPoints=V_inner; %region points
inputStruct.regionA=regionA*volumeFactor; %Desired volume for tets
inputStruct.minRegionMarker=2; %Minimum region marker

%% 
% Mesh model using tetrahedral elements using tetGen

[meshOutput]=runTetGen(inputStruct); %Run tetGen 

%%
% Visualize mesh

meshView(meshOutput);

%% 
% Access model element and patch data 
F=meshOutput.faces;
V=meshOutput.nodes;
C=meshOutput.faceMaterialID;
E=meshOutput.elements;
elementMaterialID=meshOutput.elementMaterialID;

Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;

%% Define boundary condition node sets

logicRigid= Cb==1 | Cb==3;
bcSupportList=Fb(logicRigid,:);
bcSupportList=unique(bcSupportList(:));

logicIndentor= Cb==4;
bcPrescribeList=Fb(logicIndentor,:);
bcPrescribeList=unique(bcPrescribeList(:));

%%
% Visualize boundary conditions

cFigure; hold on;
gpatch(Fb,V,'bw','none',0.5);
hp(1)=plotV(V(bcSupportList,:),'k.','lineWidth',lineWidth1,'MarkerSize',markerSize1);
hp(2)=plotV(V(bcPrescribeList,:),'r.','lineWidth',lineWidth1,'MarkerSize',markerSize1);
legend(hp,{'BC Full support','BC Prescribed displacement'});
axisGeom(gca,fontSize);
drawnow; 

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='2.5'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

%Create control structure for use by all steps
stepStruct.Control.analysis.ATTR.type='static';
stepStruct.Control.time_steps=numTimeSteps;
stepStruct.Control.step_size=1/numTimeSteps;
stepStruct.Control.time_stepper.dtmin=dtmin;
stepStruct.Control.time_stepper.dtmax=dtmax; 
stepStruct.Control.time_stepper.max_retries=max_retries;
stepStruct.Control.time_stepper.opt_iter=opt_iter;
stepStruct.Control.max_refs=max_refs;
stepStruct.Control.max_ups=max_ups;

%Add template based default settings to proposed control section
[stepStruct.Control]=structComplete(stepStruct.Control,febio_spec.Control,1); %Complement provided with default if missing

%Remove control field (part of template) since step specific control sections are used
febio_spec=rmfield(febio_spec,'Control'); 

febio_spec.Step{1}.Control=stepStruct.Control;
febio_spec.Step{1}.ATTR.id=1;
febio_spec.Step{2}.Control=stepStruct.Control;
febio_spec.Step{2}.ATTR.id=2;

%Material section
febio_spec.Material.material{1}.ATTR.type='Ogden';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c1=c1;
febio_spec.Material.material{1}.m1=m1;
febio_spec.Material.material{1}.c2=c1;
febio_spec.Material.material{1}.m2=-m1;
febio_spec.Material.material{1}.k=k;

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='tet4'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name='Tissue'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=E;

% -> NodeSets
febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcSupportList(:);

febio_spec.Geometry.NodeSet{2}.ATTR.name='bcPrescribeList';
febio_spec.Geometry.NodeSet{2}.node.ATTR.id=bcPrescribeList(:);

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.fix{1}.ATTR.bc='x';
febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{2}.ATTR.bc='y';
febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{3}.ATTR.bc='z';
febio_spec.Boundary.fix{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;

% -> Prescribe boundary conditions
%STEP 1 Up/down
febio_spec.Step{1}.Boundary.prescribe{1}.ATTR.bc='x';
febio_spec.Step{1}.Boundary.prescribe{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Step{1}.Boundary.prescribe{1}.scale.ATTR.lc=1;
febio_spec.Step{1}.Boundary.prescribe{1}.scale.VAL=1;
febio_spec.Step{1}.Boundary.prescribe{1}.relative=1;
febio_spec.Step{1}.Boundary.prescribe{1}.value=0;

febio_spec.Step{1}.Boundary.prescribe{2}.ATTR.bc='y';
febio_spec.Step{1}.Boundary.prescribe{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Step{1}.Boundary.prescribe{2}.scale.ATTR.lc=1;
febio_spec.Step{1}.Boundary.prescribe{2}.scale.VAL=1;
febio_spec.Step{1}.Boundary.prescribe{2}.relative=1;
febio_spec.Step{1}.Boundary.prescribe{2}.value=0;

febio_spec.Step{1}.Boundary.prescribe{3}.ATTR.bc='z';
febio_spec.Step{1}.Boundary.prescribe{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Step{1}.Boundary.prescribe{3}.scale.ATTR.lc=1;
febio_spec.Step{1}.Boundary.prescribe{3}.scale.VAL=1;
febio_spec.Step{1}.Boundary.prescribe{3}.relative=1;
febio_spec.Step{1}.Boundary.prescribe{3}.value=displacementMagnitude;

%STEP 3 Sideways
febio_spec.Step{2}.Boundary.prescribe{1}.ATTR.bc='x';
febio_spec.Step{2}.Boundary.prescribe{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Step{2}.Boundary.prescribe{1}.scale.ATTR.lc=2;
febio_spec.Step{2}.Boundary.prescribe{1}.scale.VAL=1;
febio_spec.Step{2}.Boundary.prescribe{1}.relative=1;
febio_spec.Step{2}.Boundary.prescribe{1}.value=displacementMagnitude;

febio_spec.Step{2}.Boundary.prescribe{2}.ATTR.bc='y';
febio_spec.Step{2}.Boundary.prescribe{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Step{2}.Boundary.prescribe{2}.scale.ATTR.lc=2;
febio_spec.Step{2}.Boundary.prescribe{2}.scale.VAL=1;
febio_spec.Step{2}.Boundary.prescribe{2}.relative=1;
febio_spec.Step{2}.Boundary.prescribe{2}.value=0;

febio_spec.Step{2}.Boundary.prescribe{3}.ATTR.bc='z';
febio_spec.Step{2}.Boundary.prescribe{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Step{2}.Boundary.prescribe{3}.scale.ATTR.lc=2;
febio_spec.Step{2}.Boundary.prescribe{3}.scale.VAL=1;
febio_spec.Step{2}.Boundary.prescribe{3}.relative=1;
febio_spec.Step{2}.Boundary.prescribe{3}.value=0;

%LoadData section
febio_spec.LoadData.loadcurve{1}.ATTR.id=1;
febio_spec.LoadData.loadcurve{1}.ATTR.type='linear';
febio_spec.LoadData.loadcurve{1}.point.VAL=[0 0; 0.25 1; 0.5 0; 0.75 -1; 1 0];

febio_spec.LoadData.loadcurve{2}.ATTR.id=2;
febio_spec.LoadData.loadcurve{2}.ATTR.type='linear';
febio_spec.LoadData.loadcurve{2}.point.VAL=[1 0; 1.25 1; 1.5 0; 1.75 -1; 2 0];

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

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_strain;
febio_spec.Output.logfile.element_data{1}.ATTR.data='E1;E2;E3';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.element_data{1}.VAL=1:size(E,1);

%% Quick viewing of the FEBio input file structure
% The |febView| function can be used to view the xml structure in a MATLAB
% figure window. 

%%
% |febView(febio_spec); %Viewing the febio file|

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
febioAnalysis.runMode='internal';%'internal' or 'external';
febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results 

if runFlag==1 %i.e. a succesful run
    
    % Importing nodal displacements from a log file
    [time_mat, N_disp_mat,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp)); %Nodal displacements    
    time_mat=[0; time_mat(:)]; %Time
    
    N_disp_mat=N_disp_mat(:,2:end,:);
    sizImport=size(N_disp_mat);
    sizImport(3)=sizImport(3)+1;
    N_disp_mat_n=zeros(sizImport);
    N_disp_mat_n(:,:,2:end)=N_disp_mat;
    N_disp_mat=N_disp_mat_n;
    
    DN_MAG=sqrt(sum(N_disp_mat.^2,2));
    DN=N_disp_mat(:,:,end);
    DN_magnitude=sqrt(sum(DN(:,3).^2,2));
    V_def=V+DN;
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
    X_DEF=V_DEF(:,1,:);
    Y_DEF=V_DEF(:,2,:);
    Z_DEF=V_DEF(:,3,:);
    [CF_def]=vertexToFaceMeasure(Fb,DN_magnitude);
    
     %%
    % Importing element data from a log file
    [~,E_data,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_strain)); %Element data
    
    %Remove nodal index column
    E_data=E_data(:,2:end,:);
    
    E_data=sqrt(0.5*( (E_data(:,1,:)-E_data(:,2,:)).^2 + (E_data(:,2,:)-E_data(:,3,:)).^2 + (E_data(:,3,:)-E_data(:,1,:)).^2 ));
        
    %Add initial state i.e. zero 
    sizImport=size(E_data); 
    sizImport(3)=sizImport(3)+1;
    E_data_mat_n=zeros(sizImport);
    E_data_mat_n(:,:,2:end)=E_data;
    E_data=E_data_mat_n;
        
    VE=patchCentre(E,V);
    logicCutElements=VE(:,2)>=0;
    
    [F_cut,CF_cut_data]=element2patch(E(logicCutElements,:),E_data(logicCutElements,:,1));
    [indBoundary]=tesBoundary(F_cut,V);
    
    CV=faceToVertexMeasure(F_cut(indBoundary,:),V,CF_cut_data(indBoundary,:));
    
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure      
    gtitle([febioFebFileNamePart,': Press play to animate']);
%     subplot(1,2,1);    
    hp1=gpatch(Fb,V_def,'kw','none',0.25); %Add graphics object to animate
    hp2=gpatch(F_cut(indBoundary,:),V_def,CV,'k',1); %Add graphics object to animate
    hp2.FaceColor='interp';
%     gpatch(Fb,V,0.5*ones(1,3),'k',0.25); %A static graphics object    
    colormap(gjet(250)); colorbar;
    caxis([0 max(E_data(:))/3]);
    axisGeom(gca,fontSize);
    axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);
    axis manual; 
    camlight headlight;   
    
    drawnow;     
    
    % Set up animation features
    animStruct.Time=time_mat; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        DN=N_disp_mat(:,:,qt); %Current displacement
        DN_magnitude=sqrt(sum(DN.^2,2)); %Current displacement magnitude
        V_def=V+DN; %Current nodal coordinates
%         [CF_def]=vertexToFaceMeasure(Fb,DN_magnitude); %Current color data to use
    
        [~,CF_cut_data]=element2patch(E(logicCutElements,:),E_data(logicCutElements,:,qt));
        CV=faceToVertexMeasure(F_cut(indBoundary,:),V,CF_cut_data(indBoundary,:));
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1 hp2 hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,V_def,CV}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
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
% Copyright (C) 2019  Kevin Mattheus Moerman
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
