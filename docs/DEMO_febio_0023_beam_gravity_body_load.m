%% DEMO_febio_0023_beam_gravity_body_load
% Below is a demonstration for:
% 
% * Building geometry for a beam with hexahedral elements
% * Defining the boundary conditions 
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement results

%% Keywords
%
% * febio_spec version 2.5
% * febio, FEBio
% * beam, gravity, body load
% * surface pressure boundary condition
% * hexahedral elements, hex8
% * beam, rectangular
% * static, solid
% * hyperelastic, Ogden
% * displacement logfile

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
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacements

%Specifying dimensions and number of elements
beamWidth=0.1; 
sampleWidth=beamWidth; %Width 
sampleThickness=5*beamWidth; %Thickness 
sampleHeight=beamWidth; %Height
pointSpacings=0.025*ones(1,3); %Desired point spacing between nodes
numElementsWidth=round(sampleWidth/pointSpacings(1)); %Number of elemens in dir 1
numElementsThickness=round(sampleThickness/pointSpacings(2)); %Number of elemens in dir 2
numElementsHeight=round(sampleHeight/pointSpacings(3)); %Number of elemens in dir 3

%Material parameter set
c1=1000; %Shear-modulus-like parameter
m1=8; %Material parameter setting degree of non-linearity
k_factor=1e2; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus

%Determine load
sampleVolume=sampleWidth.*sampleThickness.*sampleHeight;
sampleWeight=(sampleVolume.*1000);
gravityVector=[0 0 9.81.*sampleWeight];

% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
min_residual=1e-40;

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

%Define supported node set
logicFace=Cb==4; %Logic for current face set
Fr=Fb(logicFace,:); %The current face set
bcSupportList=unique(Fr(:)); %Node set part of selected face

%% 
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists. 

hf=cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,V,'kw','k',0.5);

hl(1)=plotV(V(bcSupportList,:),'k.','MarkerSize',markerSize);

legend(hl,{'BC support'});

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

%Control section
febio_spec.Control.analysis.ATTR.type='static';
febio_spec.Control.title='Cube analysis';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;
febio_spec.Control.max_refs=max_refs;
febio_spec.Control.max_ups=max_ups;
febio_spec.Control.min_residual=min_residual;

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
febio_spec.Geometry.Elements{1}.ATTR.type='hex8'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name='Beam'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=E;

% -> NodeSets
febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcSupportList(:);

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.fix{1}.ATTR.bc='x';
febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{2}.ATTR.bc='y';
febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{3}.ATTR.bc='z';
febio_spec.Boundary.fix{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;

% -> Surface load
febio_spec.Loads.body_load{1}.ATTR.type='const';
febio_spec.Loads.body_load{1}.x.VAL=gravityVector(1);
febio_spec.Loads.body_load{1}.x.ATTR.lc=1;
febio_spec.Loads.body_load{1}.y.VAL=gravityVector(2);
febio_spec.Loads.body_load{1}.y.ATTR.lc=1;
febio_spec.Loads.body_load{1}.z.VAL=gravityVector(3);
febio_spec.Loads.body_load{1}.z.ATTR.lc=1;

%Output section 
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

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
