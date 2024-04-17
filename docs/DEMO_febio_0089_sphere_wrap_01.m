%% DEMO_febio_0089_sphere_wrap_01
% Below is a demonstration for:
% 
% * Building geometry for a slab with hexahedral elements, and a
% triangulated sphere. 
% * Defining the boundary conditions 
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement results

%% Keywords
%
% * febio_spec version 4.0
% * febio, FEBio
% * indentation
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

%% Plot settings
fontSize=15;
faceAlpha1=0.8;
faceAlpha2=0.3;
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
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress
febioLogFileName_contactPressure=[febioFebFileNamePart,'_contactPressure_out.txt']; %Log file name for exporting contact pressure

%Sphere parameters
numRefineStepsSphere = 3; 
sphereRadius = 1;
f=1.001; 
layerThickness=0.01;

%Material parameter set
c1=1e-3; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
k_factor=100; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus

% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=10; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
symmetric_stiffness=0;

runMode='external';% 'internal' or 'external'

%Contact parameters
contactPenalty=10;
laugon=0;
minaug=1;
maxaug=10;
fric_coeff=0.01;

%% Creating model geometry and mesh

% Building spherical surface models
[E2,V2,~]=geoSphere(numRefineStepsSphere+1,sphereRadius*f);
E2=fliplr(E2);

[F,V,~]=geoSphere(numRefineStepsSphere,sphereRadius);
C=ones(size(F,1),1);
V_regions=[0 0 0];
V_holes=[];
regionTetVolumes=tetVolMeanEst(F,V);
stringOpt='-pq1.2AaY';

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
inputStruct.tetType='tet4'; %Set desired element type

% Mesh model using tetrahedral elements using tetGen 
[meshStruct]=runTetGen(inputStruct); %Run tetGen 

%% 
% Access mesh output structure

E1=meshStruct.elements; %The elements
V1=meshStruct.nodes; %The vertices or nodes
Fb1=meshStruct.facesBoundary; %The boundary faces
Cb1=meshStruct.boundaryMarker; %The boundary markers

%% 
% Plotting model boundary surfaces and a cut view

hFig=cFigure; 

subplot(1,2,1); hold on; 
title('Model boundary surfaces and labels','FontSize',fontSize);
gpatch(Fb1,V1,Cb1,'k',faceAlpha1); 

gpatch(E2,V2,'kw','k',faceAlpha1); 
colormap(gjet(6)); icolorbar;
axisGeom(gca,fontSize);
camlight headlight; 

hs=subplot(1,2,2); hold on; 
title('Cut view of solid mesh','FontSize',fontSize);
optionStruct.hFig=[hFig hs];
gpatch(E2,V2,'kw','k',0.5); 
meshView(meshStruct,optionStruct);
axisGeom(gca,fontSize);

drawnow;

%% Joining node sets
V=[V1;V2;]; %Combined node sets
E2=E2+size(V1,1); %Fixed element indices

%%
% Plotting joined geometry
cFigure;
title('Joined node sets','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
gpatch(Fb1,V,Cb1,'k',faceAlpha1); 
gpatch(E2,V,'kw','k',faceAlpha1);
colormap(gjet(6)); icolorbar; 
axisGeom(gca,fontSize);
camlight headlight;
drawnow;

%% Define contact surfaces

% The rigid primary surface of the sphere
F_contact_secondary=E2;

% The deformable secondary surface of the slab
F_contact_primary=fliplr(Fb1);

% Plotting surface models
cFigure; hold on;
title('Contact sets and normal directions','FontSize',fontSize);

gpatch(Fb1,V,'kw','none',faceAlpha2); 
hl(1)=gpatch(F_contact_secondary,V,'g','k',0.5); 
patchNormPlot(F_contact_secondary,V);
hl(2)=gpatch(F_contact_primary,V,'b','k',0.5);
patchNormPlot(F_contact_primary,V);

legend(hl,{'Secondary','Primary'});

axisGeom(gca,fontSize);
camlight headlight;
drawnow;

%% Define boundary conditions

V_final = V; 
V_final(:,3)= V_final(:,3)/2; 
U = V_final-V; 

cFigure; hold on;
gpatch(E2,V+U,'kw','k',faceAlpha2); 

% quiverVec(V,U)
axisGeom(gca,fontSize);
camlight headlight;
drawnow;


% 
% %Supported nodes
% logicRigid=Cb1==5;
% Fr=Fb1(logicRigid,:);
% bcSupportList=unique(Fr(:));
% 
% %%
% % Visualize BC's
% hf=cFigure;
% title('Boundary conditions model','FontSize',fontSize);
% xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
% hold on;
% 
% gpatch(Fb1,V,'kw','none',faceAlpha2); 
% hl2(1)=gpatch(E2,V,'kw','k',1); 
% 
% hl2(2)=plotV(V(bcSupportList,:),'k.','MarkerSize',markerSize);
% 
% legend(hl2,{'Rigid body sphere','BC support'});
% 
% axisGeom(gca,fontSize);
% camlight headlight;
% drawnow;

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
% febio_spec.Material.material{1}.c2=c1;
% febio_spec.Material.material{1}.m2=-m1;
febio_spec.Material.material{1}.k=k;

materialName2='Material2';
febio_spec.Material.material{2}.ATTR.name=materialName2;
febio_spec.Material.material{2}.ATTR.type='Ogden';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.c1=c1;
febio_spec.Material.material{2}.m1=m1;
% febio_spec.Material.material{2}.c2=c1;
% febio_spec.Material.material{2}.m2=-m1;
febio_spec.Material.material{2}.k=k;

%Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='tet4'; %Element type 
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E1,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E1; %The element matrix

partName2='Part2';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='tri3'; %Element type 
febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=E2; %The element matrix

% -> NodeSets
nodeSetName1='bcSupportList';
bcPrescribeList=unique(E2(:));
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcPrescribeList);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat=materialName1;

febio_spec.MeshDomains.ShellDomain{1}.ATTR.name=partName2;
febio_spec.MeshDomains.ShellDomain{1}.ATTR.mat=materialName2;
febio_spec.MeshDomains.ShellDomain{1}.shell_thickness=layerThickness; 

% -> Surfaces
surfaceName1='contactSurface1';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.tri3.ATTR.id=(1:1:size(F_contact_primary,1))';
febio_spec.Mesh.Surface{1}.tri3.VAL=F_contact_primary;

surfaceName2='contactSurface2';
febio_spec.Mesh.Surface{2}.ATTR.name=surfaceName2;
febio_spec.Mesh.Surface{2}.tri3.ATTR.id=(1:1:size(F_contact_secondary,1))';
febio_spec.Mesh.Surface{2}.tri3.VAL=F_contact_secondary;

% -> Surface pairs
contactPairName='Contact1';
febio_spec.Mesh.SurfacePair{1}.ATTR.name=contactPairName;
febio_spec.Mesh.SurfacePair{1}.primary=surfaceName1;
febio_spec.Mesh.SurfacePair{1}.secondary=surfaceName2;

loadDataName1='nodal_displacements';
febio_spec.MeshData.NodeData.ATTR.name=loadDataName1;
febio_spec.MeshData.NodeData.ATTR.node_set=nodeSetName1;
febio_spec.MeshData.NodeData.ATTR.data_type='vec3';
febio_spec.MeshData.NodeData.node.ATTR.lid=(1:1:numel(bcPrescribeList))';
febio_spec.MeshData.NodeData.node.VAL=U(bcPrescribeList,:);

% Boundary condition section 

febio_spec.Boundary.bc{1}.ATTR.name='bcPrescribeList01_01';
febio_spec.Boundary.bc{1}.ATTR.type='prescribed displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.dof='x';
febio_spec.Boundary.bc{1}.value.ATTR.lc=1;
febio_spec.Boundary.bc{1}.value.VAL=0.5;
febio_spec.Boundary.bc{1}.relative=1;

febio_spec.Boundary.bc{2}.ATTR.name='bcPrescribeList01_02';
febio_spec.Boundary.bc{2}.ATTR.type='prescribed displacement';
febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{2}.dof='y';
febio_spec.Boundary.bc{2}.value.ATTR.lc=1;
febio_spec.Boundary.bc{2}.value.VAL=0;
febio_spec.Boundary.bc{2}.relative=1;

febio_spec.Boundary.bc{3}.ATTR.name='bcPrescribeList01_03';
febio_spec.Boundary.bc{3}.ATTR.type='prescribed displacement';
febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{3}.dof='z';
febio_spec.Boundary.bc{3}.value.ATTR.lc=1;
febio_spec.Boundary.bc{3}.value.VAL=0;
febio_spec.Boundary.bc{3}.relative=1;



% % -> Fix boundary conditions
% febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_xyz';
% febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
% febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
% febio_spec.Boundary.bc{1}.x_dof=1;
% febio_spec.Boundary.bc{1}.y_dof=1;
% febio_spec.Boundary.bc{1}.z_dof=1;

%Contact section
febio_spec.Contact.contact{1}.ATTR.type='sliding-elastic';
febio_spec.Contact.contact{1}.ATTR.surface_pair=contactPairName;
febio_spec.Contact.contact{1}.two_pass=0;
febio_spec.Contact.contact{1}.laugon=laugon;
febio_spec.Contact.contact{1}.tolerance=0.2;
febio_spec.Contact.contact{1}.gaptol=0;
febio_spec.Contact.contact{1}.minaug=minaug;
febio_spec.Contact.contact{1}.maxaug=maxaug;
febio_spec.Contact.contact{1}.search_tol=0.01;
febio_spec.Contact.contact{1}.search_radius=0.1*sqrt(sum((max(V,[],1)-min(V,[],1)).^2,2));
febio_spec.Contact.contact{1}.symmetric_stiffness=0;
febio_spec.Contact.contact{1}.auto_penalty=1;
febio_spec.Contact.contact{1}.update_penalty=1;
febio_spec.Contact.contact{1}.penalty=contactPenalty;
febio_spec.Contact.contact{1}.fric_coeff=fric_coeff;

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

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='s3';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.face_data{1}.ATTR.file=febioLogFileName_contactPressure;
febio_spec.Output.logfile.face_data{1}.ATTR.data='contact pressure';
febio_spec.Output.logfile.face_data{1}.ATTR.surface=surfaceName1;
febio_spec.Output.logfile.face_data{1}.ATTR.delim=',';

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
febioAnalysis.runMode=runMode;

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results 

if runFlag==1 %i.e. a succesful run
    
    %% 
    % Importing nodal displacements from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),0,1);
    
    %Access data
    N_disp_mat=dataStruct.data; %Displacement
    timeVec=dataStruct.time; %Time
    
    %Create deformed coordinate set
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
            
    %%
    % Importing element stress from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),0,1);     
    
    %Access data
    E_stress_mat=dataStruct.data;
    
                
    %%
    % Importing contact pressure from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_contactPressure),0,1);
    
    %Access data
    F_primary_contact_pressure_mat=dataStruct.data;
    
    %%
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations
    
    [CV]=faceToVertexMeasure(E1,V,E_stress_mat(:,:,end));
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('$\sigma_{3}$ [MPa]','Interpreter','Latex')
    hp=gpatch(Fb1,V_DEF(:,:,end),CV,'k',1); %Add graphics object to animate
    hp.Marker='.';
    hp.MarkerSize=markerSize2;
    hp.FaceColor='interp';
    
    hp2=gpatch(E2,V_DEF(:,:,end),'w','none',0.5); %Add graphics object to animate
    
    axisGeom(gca,fontSize);
    colormap(flipud(gjet(250))); colorbar;
    caxis([min(E_stress_mat(:)) max(E_stress_mat(:))]);
    axis(axisLim(V_DEF)); %Set axis limits statically
    camlight headlight;
    
    % Set up animation features
    animStruct.Time=timeVec; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments
        
        [CV]=faceToVertexMeasure(E1,V,E_stress_mat(:,:,qt));
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),CV,V_DEF(:,:,qt)}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct); %Initiate animation feature
    drawnow;
    
        %%
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations

    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('Contact pressure [MPa]','Interpreter','Latex')
    gpatch(Fb1,V_DEF(:,:,end),'w','none',0.1);
    hp=gpatch(F_contact_primary,V_DEF(:,:,end),F_primary_contact_pressure_mat(:,:,1),'k',1); %Add graphics object to animate        
        
    hp2=gpatch(E2,V_DEF(:,:,end),'w','none',0.25); %Add graphics object to animate
    
    axisGeom(gca,fontSize);
    colormap(gjet(250)); colorbar;
    caxis([min(F_primary_contact_pressure_mat(:)) max(F_primary_contact_pressure_mat(:))]);
    axis(axisLim(V_DEF)); %Set axis limits statically
    camlight headlight;
    
    % Set up animation features
    animStruct.Time=timeVec; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),F_primary_contact_pressure_mat(:,:,qt),V_DEF(:,:,qt)}; %Property values for to set in order to animate
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
