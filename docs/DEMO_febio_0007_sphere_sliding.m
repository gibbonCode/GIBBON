%% DEMO_febio_0007_sphere_sliding
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
% * febio_spec version 2.5
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
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress

%Specifying dimensions and number of elements for slab
sampleHeight=5; %Height
sampleWidth=sampleHeight*4; %Width 
sampleThickness=sampleHeight*1; %Thickness 
pointSpacings=1*ones(1,3); %Desired point spacing between nodes
numElementsWidth=round(sampleWidth/pointSpacings(1)); %Number of elemens in dir 1
numElementsThickness=round(sampleThickness/pointSpacings(2)); %Number of elemens in dir 2
numElementsHeight=round(sampleHeight/pointSpacings(3)); %Number of elemens in dir 3

%Sphere parameters
numRefineStepsSphere=2; 
sphereRadius=sampleHeight/2;

%Define applied displacement
sphereIndentationDisplacement=sphereRadius; 
sphereSlideDisplacement=sampleWidth-(sphereRadius*2); 

%Material parameter set
c1=1e-3; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
k_factor=100; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus

% FEA control settings
numTimeSteps=20; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=15; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
runMode='external';%'internal';

%Contact parameters
contactInitialOffset=0.1;
contactPenalty=10;
contactAlg=2;
switch contactAlg
    case 1
        contactType='sticky';
    case 2
        contactType='facet-to-facet sliding'; 
    case 3
        contactType='sliding_with_gaps';
    case 4
        contactType='sliding2';
end

%% Creating model geometry and mesh
% A box is created with tri-linear hexahedral (hex8) elements using the
% |hexMeshBox| function. The function offers the boundary faces with
% seperate labels for the top, bottom, left, right, front, and back sides.
% As such these can be used to define boundary conditions on the exterior. 

% Create a box with hexahedral elements
beamDimensions=[sampleWidth sampleThickness sampleHeight]; %Dimensions
beamElementNumbers=[numElementsWidth numElementsThickness numElementsHeight]; %Number of elements
outputStructType=2; %A structure compatible with mesh view
[meshStruct]=hexMeshBox(beamDimensions,beamElementNumbers,outputStructType);

%Access elements, nodes, and faces from the structure
E1=meshStruct.elements; %The elements 
V1=meshStruct.nodes; %The nodes (vertices)
Fb1=meshStruct.facesBoundary; %The boundary faces
Cb1=meshStruct.boundaryMarker; %The "colors" or labels for the boundary faces
elementMaterialIndices=ones(size(E1,1),1); %Element material indices

%% Creating triangulated sphere surface model

[E2,V2,~]=geoSphere(numRefineStepsSphere,sphereRadius); 

%Offset indentor
minV2=min(V2,[],1);
minV1=min(V1,[],1);

V2(:,1)=V2(:,1)-minV2(1)+minV1(1);

V2(:,3)=V2(:,3)-minV2(3)+(sampleHeight/2)+contactInitialOffset;

center_of_mass=mean(V2,1);

%% 
% Plotting model boundary surfaces and a cut view

hFig=cFigure; 

subplot(1,2,1); hold on; 
title('Model boundary surfaces and labels','FontSize',fontSize);
gpatch(Fb1,V1,Cb1,'k',faceAlpha1); 
gpatch(E2,V2,'kw','k',faceAlpha1); 
colormap(gjet(6)); icolorbar;
axisGeom(gca,fontSize);

hs=subplot(1,2,2); hold on; 
title('Cut view of solid mesh','FontSize',fontSize);
optionStruct.hFig=[hFig hs];
gpatch(E2,V2,'kw','k',1); 
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

% The rigid master surface of the sphere
F_contact_master=E2;

% The deformable slave surface of the slab
logicContactSurf1=Cb1==6;
F_contact_slave=Fb1(logicContactSurf1,:);

% Plotting surface models
cFigure; hold on;
title('Contact sets and normal directions','FontSize',fontSize);

gpatch(Fb1,V,'kw','none',faceAlpha2); 
hl(1)=gpatch(F_contact_master,V,'g','k',1); 
patchNormPlot(F_contact_master,V);
hl(2)=gpatch(F_contact_slave,V,'b','k',1);
patchNormPlot(F_contact_slave,V);

legend(hl,{'Master','Slave'});

axisGeom(gca,fontSize);
camlight headlight;
drawnow;

%% Define boundary conditions

%Supported nodes
logicRigid=Cb1==5;
Fr=Fb1(logicRigid,:);
bcSupportList=unique(Fr(:));

%Prescribed displacement nodes
bcPrescribeList=unique(E2(:));
bcPrescribeMagnitudes=[sphereSlideDisplacement 0 -(sphereIndentationDisplacement+contactInitialOffset)];

%Visualize BC's
hf=cFigure;
title('Boundary conditions model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb1,V,'kw','none',faceAlpha2); 
hl2(1)=gpatch(E2,V,'kw','k',1); 

hl2(2)=plotV(V(bcSupportList,:),'k.','MarkerSize',markerSize);

legend(hl2,{'Rigid body sphere','BC support'});

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

%Step specific control section
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

febio_spec.Material.material{2}.ATTR.type='rigid body';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.density=1;
febio_spec.Material.material{2}.center_of_mass=center_of_mass;

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='hex8'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name='Slab'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(E1,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=E1;

febio_spec.Geometry.Elements{2}.ATTR.type='tri3'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat=2; %material index for this set 
febio_spec.Geometry.Elements{2}.ATTR.name='Sphere'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL=E2;

% -> NodeSets
febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcSupportList(:);

% -> Surfaces
febio_spec.Geometry.Surface{1}.ATTR.name='contact_master';
febio_spec.Geometry.Surface{1}.tri3.ATTR.lid=(1:1:size(F_contact_master,1))';
febio_spec.Geometry.Surface{1}.tri3.VAL=F_contact_master;

febio_spec.Geometry.Surface{2}.ATTR.name='contact_slave';
febio_spec.Geometry.Surface{2}.quad4.ATTR.lid=(1:1:size(F_contact_slave,1))';
febio_spec.Geometry.Surface{2}.quad4.VAL=F_contact_slave;

% -> Surface pairs
febio_spec.Geometry.SurfacePair{1}.ATTR.name='Contact1';
febio_spec.Geometry.SurfacePair{1}.master.ATTR.surface=febio_spec.Geometry.Surface{1}.ATTR.name;
febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.fix{1}.ATTR.bc='x';
febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{2}.ATTR.bc='y';
febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{3}.ATTR.bc='z';
febio_spec.Boundary.fix{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;

% -> Prescribed boundary conditions on the rigid body
febio_spec.Step{1}.Boundary.rigid_body{1}.ATTR.mat=2;
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='x';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='Rx';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Ry';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Rz';
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.ATTR.bc='z';
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.ATTR.lc=1;
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.VAL=bcPrescribeMagnitudes(3);

febio_spec.Step{2}.Boundary.rigid_body{1}.ATTR.mat=2;
febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='z';
febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='Rx';
febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Ry';
febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Rz';
febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.bc='x';
febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.lc=2;
febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.VAL=bcPrescribeMagnitudes(1);

%Contact section
switch contactType
    case 'sticky'
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='sticky';
        febio_spec.Contact.contact{1}.penalty=contactPenalty;
        febio_spec.Contact.contact{1}.laugon=0;
        febio_spec.Contact.contact{1}.tolerance=0.1;
        febio_spec.Contact.contact{1}.minaug=0;
        febio_spec.Contact.contact{1}.maxaug=10;
        febio_spec.Contact.contact{1}.snap_tol=0;
        febio_spec.Contact.contact{1}.max_traction=0;
        febio_spec.Contact.contact{1}.search_tolerance=0.1;
    case 'facet-to-facet sliding'
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='facet-to-facet sliding';
        febio_spec.Contact.contact{1}.penalty=contactPenalty;
        febio_spec.Contact.contact{1}.auto_penalty=1;
        febio_spec.Contact.contact{1}.two_pass=0;
        febio_spec.Contact.contact{1}.laugon=0;
        febio_spec.Contact.contact{1}.tolerance=0.1;
        febio_spec.Contact.contact{1}.gaptol=0;
        febio_spec.Contact.contact{1}.minaug=0;
        febio_spec.Contact.contact{1}.maxaug=10;
        febio_spec.Contact.contact{1}.search_tol=0.01;
        febio_spec.Contact.contact{1}.search_radius=mean(pointSpacings)/2;
    case 'sliding_with_gaps'
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='sliding_with_gaps';
        febio_spec.Contact.contact{1}.penalty=contactPenalty;
        febio_spec.Contact.contact{1}.auto_penalty=1;
        febio_spec.Contact.contact{1}.two_pass=0;
        febio_spec.Contact.contact{1}.laugon=0;
        febio_spec.Contact.contact{1}.tolerance=0.1;
        febio_spec.Contact.contact{1}.gaptol=0;
        febio_spec.Contact.contact{1}.minaug=0;
        febio_spec.Contact.contact{1}.maxaug=10;
        febio_spec.Contact.contact{1}.fric_coeff=0;
        febio_spec.Contact.contact{1}.fric_penalty=0;
        febio_spec.Contact.contact{1}.ktmult=1;
        febio_spec.Contact.contact{1}.seg_up=0;
        febio_spec.Contact.contact{1}.search_tol=0.01;
    case 'sliding2'
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='sliding2';
        febio_spec.Contact.contact{1}.penalty=contactPenalty;
        febio_spec.Contact.contact{1}.auto_penalty=1;
        febio_spec.Contact.contact{1}.two_pass=0;
        febio_spec.Contact.contact{1}.laugon=0;
        febio_spec.Contact.contact{1}.tolerance=0.1;
        febio_spec.Contact.contact{1}.gaptol=0;
        febio_spec.Contact.contact{1}.symmetric_stiffness=0;
        febio_spec.Contact.contact{1}.search_tol=0.01;
        febio_spec.Contact.contact{1}.search_radius=mean(pointSpacings)/2;
end

% LoadData section
febio_spec.LoadData.loadcurve{1}.ATTR.id=1;
febio_spec.LoadData.loadcurve{1}.ATTR.type='linear';
febio_spec.LoadData.loadcurve{1}.point.VAL=[0 0; 1 1; 2 1];

febio_spec.LoadData.loadcurve{2}.ATTR.id=2;
febio_spec.LoadData.loadcurve{2}.ATTR.type='linear';
febio_spec.LoadData.loadcurve{2}.point.VAL=[0 0; 1 0; 2 1];

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
febio_spec.Output.logfile.element_data{1}.ATTR.data='s3';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.element_data{1}.VAL=1:size(E1,1);

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
febioAnalysis.runMode=runMode;
febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
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
    % Importing element stress from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),1,1);     
    
    %Access data
    E_stress_mat=dataStruct.data;
    E_stress_mat(isnan(E_stress_mat))=0;
    
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
