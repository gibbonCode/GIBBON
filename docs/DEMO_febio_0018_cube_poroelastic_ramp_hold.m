%% DEMO_febio_0018_cube_poroelastic_ramp_hold
% Below is a demonstration for:
% 
% * Building geometry for a cube with hexahedral elements
% * Defining the boundary conditions 
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement and stress results

%% Keywords
%
% * febio_spec version 2.5
% * febio, FEBio
% * uniaxial loading
% * compression, tension, compressive, tensile
% * displacement control, displacement boundary condition
% * hexahedral elements, hex8
% * cube, box, rectangular
% * static, solid
% * hyperelastic, Ogden
% * viscoelastic
% * poroelastic
% * biphasic
% * uncoupled, coupled
% * ramp hold
% * displacement logfile
% * stress logfile

%%

clear; close all; clc;

%% Plot settings
fontSize=20;
faceAlpha1=0.8;
markerSize=40;
markerSize2=20;
lineWidth=3;

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress

%Specifying dimensions and number of elements
meshType='hex8'; %hex8 or tet4
unitSystem=2; %1=m, 2=mm
switch unitSystem
    case 1
        min_residual=1e-40;
        cubeSize=10e-3;
        sampleWidth=cubeSize; %Width
        sampleThickness=cubeSize; %Thickness
        sampleHeight=cubeSize; %Height
        pointSpacings=2e-3*ones(1,3); %Desired point spacing between nodes
        numElementsWidth=round(sampleWidth/pointSpacings(1)); %Number of elemens in dir 1
        numElementsThickness=round(sampleThickness/pointSpacings(2)); %Number of elemens in dir 2
        numElementsHeight=round(sampleHeight/pointSpacings(3)); %Number of elemens in dir 3
        
        %Define applied displacement
        appliedStrain=0.3; %Linear strain (Only used to compute applied stretch)
        loadingOption='compression'; % or 'tension'
        switch loadingOption
            case 'compression'
                stretchLoad=1-appliedStrain; %The applied stretch for uniaxial loading
            case 'tension'
                stretchLoad=1+appliedStrain; %The applied stretch for uniaxial loading
        end
        displacementMagnitude=(stretchLoad*sampleHeight)-sampleHeight; %The displacement magnitude
        
        %Material parameter set
        
        %Hyperelastic parameters
        c1=1000; %ogden c1
        m1=6; %ogden m1
        k_factor=1; %Bulk like modulus factor
        k=c1*k_factor; %The bulk like modulus
        
        d=1000; %Density
        
        %Constant Isotropic Permeability parameters
        phi0=0.5; %Solid volume fraction in reference configuration
        permHydro=7.41e-11; %hydraulic permeability
    case 2
        min_residual=1e-40;
        cubeSize=10;
        sampleWidth=cubeSize; %Width
        sampleThickness=cubeSize; %Thickness
        sampleHeight=cubeSize; %Height
        pointSpacings=2*ones(1,3); %Desired point spacing between nodes
        numElementsWidth=round(sampleWidth/pointSpacings(1)); %Number of elemens in dir 1
        numElementsThickness=round(sampleThickness/pointSpacings(2)); %Number of elemens in dir 2
        numElementsHeight=round(sampleHeight/pointSpacings(3)); %Number of elemens in dir 3
        
        %Define applied displacement
        appliedStrain=0.3; %Linear strain (Only used to compute applied stretch)
        loadingOption='compression'; % or 'tension'
        switch loadingOption
            case 'compression'
                stretchLoad=1-appliedStrain; %The applied stretch for uniaxial loading
            case 'tension'
                stretchLoad=1+appliedStrain; %The applied stretch for uniaxial loading
        end
        displacementMagnitude=(stretchLoad*sampleHeight)-sampleHeight; %The displacement magnitude
        
        %Material parameter set
        
        %Hyperelastic parameters
        c1=1e-3; %ogden c1
        m1=6; %ogden m1
        k_factor=1; %Bulk like modulus factor
        k=c1*k_factor; %The bulk like modulus
        
        d=1e-9; %Density
        
        %Constant Isotropic Permeability parameters
        phi0=0.5; %Solid volume fraction in reference configuration
        permHydro=7.41e1; %hydraulic permeability
end

% FEA control settings
analysisType='transient'; 
febioModule='biphasic';

t_load=0.1; %Time from start to max load
t_step_ini1=t_load/50; %Initial desired step size
numTimeSteps1=round(t_load/t_step_ini1); %Number of time steps desired
t_step1=t_load/numTimeSteps1; %Step size
dtmin1=t_step1/100; %Smallest allowed step size
dtmax1=t_step1; %Largest allowed step size

t_hold=5;
t_step_ini2=t_step_ini1; %Initial desired step size
numTimeSteps2=round(t_hold/t_step_ini2); %Number of time steps desired
t_step2=t_hold/numTimeSteps2; %Step size
dtmin2=t_step2/100; %Smallest allowed step size
dtmax2=0.25; %Largest allowed step size

max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=5; %Maximum number of retires

%% Creating model geometry and mesh
% A box is created with tri-linear hexahedral (hex8) elements using the
% |hexMeshBox| function. The function offers the boundary faces with
% seperate labels for the top, bottom, left, right, front, and back sides.
% As such these can be used to define boundary conditions on the exterior. 

% Create a box with hexahedral elements
cubeDimensions=[sampleWidth sampleThickness sampleHeight]; %Dimensions
cubeElementNumbers=[numElementsWidth numElementsThickness numElementsHeight]; %Number of elements

switch meshType
    case 'hex8' %hex8 structured
        [meshStruct]=hexMeshBox(cubeDimensions,cubeElementNumbers,2);
    case 'tet4' %tet4 unstructured
        [meshStruct]=tetMeshBox(cubeDimensions,mean(pointSpacings));
end

%Access elements, nodes, and faces from the structure
E=meshStruct.elements; %The elements 
V=meshStruct.nodes; %The nodes (vertices)
Fb=meshStruct.facesBoundary; %The boundary faces
Cb=meshStruct.boundaryMarker; %The "colors" or labels for the boundary faces
elementMaterialIndices=ones(size(E,1),1); %Element material indices

%% 
% Plotting model boundary surfaces and a cut view

hFig=cFigure; 

subplot(1,2,1); hold on; 
title('Model boundary surfaces and labels','FontSize',fontSize);
gpatch(Fb,V,Cb,'k',faceAlpha1); 
colormap(gjet(6)); icolorbar;
axisGeom(gca,fontSize);

hs=subplot(1,2,2); hold on; 
title('Cut view of solid mesh','FontSize',fontSize);
optionStruct.hFig=[hFig hs];
meshView(meshStruct,optionStruct);
axisGeom(gca,fontSize);

drawnow;

%% Defining the boundary conditions
% The visualization of the model boundary shows colors for each side of the
% cube. These labels can be used to define boundary conditions. 

%Define supported node sets
logicFace=Cb==1; %Logic for current face set
Fr=Fb(logicFace,:); %The current face set
bcSupportList_X=unique(Fr(:)); %Node set part of selected face

logicFace=Cb==3; %Logic for current face set
Fr=Fb(logicFace,:); %The current face set
bcSupportList_Y=unique(Fr(:)); %Node set part of selected face

logicFace=Cb==5; %Logic for current face set
Fr=Fb(logicFace,:); %The current face set
bcSupportList_Z=unique(Fr(:)); %Node set part of selected face

%Prescribed displacement nodes
logicPrescribe=Cb==6; %Logic for current face set
Fr=Fb(logicPrescribe,:); %The current face set
bcPrescribeList=unique(Fr(:)); %Node set part of selected face

%% 
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists. 

hf=cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,V,'kw','k',0.5);

hl(1)=plotV(V(bcSupportList_X,:),'r.','MarkerSize',markerSize);
hl(2)=plotV(V(bcSupportList_Y,:),'g.','MarkerSize',markerSize);
hl(3)=plotV(V(bcSupportList_Z,:),'b.','MarkerSize',markerSize);
hl(4)=plotV(V(bcPrescribeList,:),'k.','MarkerSize',markerSize);

legend(hl,{'BC x support','BC y support','BC z support','BC z prescribe'});

axisGeom(gca,fontSize);
camlight headlight; 
drawnow; 

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='2.5'; 

%Module section
febio_spec.Module.ATTR.type=febioModule; 

%Get control section from template
stepStruct.Control=febio_spec.Control;
stepStruct.Control.symmetric_stiffness=0; %Recommended for biphasic analysis

%Remove control field (part of template) since step specific control sections are used
febio_spec=rmfield(febio_spec,'Control'); 

%Control sections for each step
febio_spec.Step{1}.ATTR.id=1;
febio_spec.Step{1}.Control=stepStruct.Control;
febio_spec.Step{1}.Control.analysis.ATTR.type=analysisType;
febio_spec.Step{1}.Control.time_steps=numTimeSteps1;
febio_spec.Step{1}.Control.step_size=t_step1;
febio_spec.Step{1}.Control.time_stepper.dtmin=dtmin1;
febio_spec.Step{1}.Control.time_stepper.dtmax=dtmax1; 
febio_spec.Step{1}.Control.time_stepper.max_retries=max_retries;
febio_spec.Step{1}.Control.time_stepper.opt_iter=opt_iter;
febio_spec.Step{1}.Control.max_refs=max_refs;
febio_spec.Step{1}.Control.max_ups=max_ups;
febio_spec.Step{1}.Control.min_residual=min_residual;

febio_spec.Step{2}.ATTR.id=2;
febio_spec.Step{2}.Control=stepStruct.Control;
febio_spec.Step{2}.Control.analysis.ATTR.type=analysisType;
febio_spec.Step{2}.Control.time_steps=numTimeSteps2;
febio_spec.Step{2}.Control.step_size=t_step2;
febio_spec.Step{2}.Control.time_stepper.dtmin=dtmin2;
febio_spec.Step{2}.Control.time_stepper.dtmax=dtmax2; 
febio_spec.Step{2}.Control.time_stepper.max_retries=max_retries;
febio_spec.Step{2}.Control.time_stepper.opt_iter=opt_iter;
febio_spec.Step{2}.Control.max_refs=max_refs;
febio_spec.Step{2}.Control.max_ups=max_ups;
febio_spec.Step{2}.Control.min_residual=min_residual;

%Material section

%Viscous part
febio_spec.Material.material{1}.ATTR.type='biphasic';
febio_spec.Material.material{1}.ATTR.name='Block_material';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.phi0=phi0;
febio_spec.Material.material{1}.permeability.ATTR.type='perm-const-iso';
febio_spec.Material.material{1}.permeability.ATTR.name='permeability';
febio_spec.Material.material{1}.permeability.perm=permHydro;
febio_spec.Material.material{1}.fluid_density=d;

%Solid part
febio_spec.Material.material{1}.solid{1}.ATTR.type='Ogden unconstrained';
febio_spec.Material.material{1}.solid{1}.c1=c1;
febio_spec.Material.material{1}.solid{1}.m1=m1;
% febio_spec.Material.material{1}.solid{1}.c2=c1;
% febio_spec.Material.material{1}.solid{1}.m2=-m1;
febio_spec.Material.material{1}.solid{1}.cp=k;
febio_spec.Material.material{1}.solid{1}.density=d;

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type=meshType; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name='Cube'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=E;

% -> NodeSets
febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList_X';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcSupportList_X(:);

febio_spec.Geometry.NodeSet{2}.ATTR.name='bcSupportList_Y';
febio_spec.Geometry.NodeSet{2}.node.ATTR.id=bcSupportList_Y(:);

febio_spec.Geometry.NodeSet{3}.ATTR.name='bcSupportList_Z';
febio_spec.Geometry.NodeSet{3}.node.ATTR.id=bcSupportList_Z(:);

febio_spec.Geometry.NodeSet{4}.ATTR.name='bcPrescribeList';
febio_spec.Geometry.NodeSet{4}.node.ATTR.id=bcPrescribeList(:);

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.fix{1}.ATTR.bc='x';
febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{2}.ATTR.bc='y';
febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Boundary.fix{3}.ATTR.bc='z';
febio_spec.Boundary.fix{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{3}.ATTR.name;

% -> Prescribe boundary conditions
febio_spec.Boundary.prescribe{1}.ATTR.bc='z';
febio_spec.Boundary.prescribe{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{4}.ATTR.name;
febio_spec.Boundary.prescribe{1}.scale.ATTR.lc=1;
febio_spec.Boundary.prescribe{1}.scale.VAL=1;
febio_spec.Boundary.prescribe{1}.relative=1;
febio_spec.Boundary.prescribe{1}.value=displacementMagnitude;

febio_spec.Boundary.prescribe{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{4}.ATTR.name;
febio_spec.Boundary.prescribe{2}.ATTR.bc='p';
febio_spec.Boundary.prescribe{2}.scale.ATTR.lc=1;
febio_spec.Boundary.prescribe{2}.scale.VAL=0;

%LoadData section 
% -> Load curves
febio_spec.LoadData.loadcurve{1}.ATTR.id=1;
febio_spec.LoadData.loadcurve{1}.ATTR.type='linear';
febio_spec.LoadData.loadcurve{1}.point.VAL=[0 0;t_load 1;(t_load+t_hold) 1];

%Output section 
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='sz';
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
febioAnalysis.runMode='external';%'internal';
febioAnalysis.t_check=0.1; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results 


if runFlag==1 %i.e. a succesful run
    
    %% 
    % Importing nodal displacements from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
    
    %Access data
    N_disp_mat=dataStruct.data; %Displacement
    timeVec=dataStruct.time; %Time
    
    %Create deformed coordinate set
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
               
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    DN_magnitude=sqrt(sum(N_disp_mat(:,:,end).^2,2)); %Current displacement magnitude
        
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('Displacement magnitude [mm]','Interpreter','Latex')
    hp=gpatch(Fb,V_DEF(:,:,end),DN_magnitude,'k',1); %Add graphics object to animate
    hp.Marker='.';
    hp.MarkerSize=markerSize2;
    hp.FaceColor='interp';
    gpatch(Fb,V,0.5*ones(1,3),'k',0.25); %A static graphics object
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]);    
    axis(axisLim(V_DEF)); %Set axis limits statically    
    camlight headlight;        
        
    % Set up animation features
    animStruct.Time=timeVec; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        DN_magnitude=sqrt(sum(N_disp_mat(:,:,qt).^2,2)); %Current displacement magnitude
                
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),DN_magnitude}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;
            
    %%
    % Importing element stress from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),1,1);
    
    %Access data
    E_stress_mat=dataStruct.data;
    
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    [CV]=faceToVertexMeasure(E,V,E_stress_mat(:,:,end));
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('$\sigma_{zz}$ [MPa]','Interpreter','Latex')
    hp=gpatch(Fb,V_DEF(:,:,end),CV,'k',1); %Add graphics object to animate
    hp.Marker='.';
    hp.MarkerSize=markerSize2;
    hp.FaceColor='interp';
    gpatch(Fb,V,0.5*ones(1,3),'k',0.25); %A static graphics object
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([min(E_stress_mat(:)) max(E_stress_mat(:))]);    
    axis(axisLim(V_DEF)); %Set axis limits statically    
    camlight headlight;        
        
    % Set up animation features
    animStruct.Time=timeVec; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        
        [CV]=faceToVertexMeasure(E,V,E_stress_mat(:,:,qt));
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;
    
    %% 
    % Calculate metrics to visualize time-stress curve
    
    stress_cauchy_sim=mean(squeeze(E_stress_mat),1)';
    
    %%    
    % Visualize stress-stretch curve
    
    cFigure; hold on;    
    title('Uniaxial stress-time curve','FontSize',fontSize);
    xlabel('Time [s]','FontSize',fontSize,'Interpreter','Latex'); 
    ylabel('$\sigma_{zz}$ [MPa]','FontSize',fontSize,'Interpreter','Latex'); 
    
    plot(timeVec(:),stress_cauchy_sim(:),'r-','lineWidth',lineWidth);
    
    view(2); axis tight;  grid on; axis square; box on; 
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
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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
