%% DEMO_febio_0015_cube_fibers_transiso
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
loadingOption='tension'; % or 'compression'
switch loadingOption
    case 'compression'
        stretchLoad=1-appliedStrain; %The applied stretch for uniaxial loading
    case 'tension'
        stretchLoad=1+appliedStrain; %The applied stretch for uniaxial loading
end
displacementMagnitude=(stretchLoad*sampleHeight)-sampleHeight; %The displacement magnitude

%Material parameter set
k_factor=100; %Bulk modulus factor    

fiberType=1; 
switch fiberType
    case 1
        c1=0.2; %Shear-modulus-like parameter
        m1=2; %Material parameter setting degree of non-linearity of ground matrix
        ksi=c1/8;
        alphaPar=2;
        beta=2;          
        k=0.5.*(c1+ksi)*k_factor; %Bulk modulus
    case 2
        Q=0.5;
        
        c=6.115e-04;
        c1=Q*c;
        m1=2.007;
        ksi_p=9.0590e-05;
        beta=3.294*ones(1,3);
        f_transiso=235.1;
        k=(2*c1+ksi_p)*k_factor;
        ksi=[ksi_p ksi_p f_transiso*ksi_p];
end
alphaFib=(45/180)*pi;

% FEA control settings
numTimeSteps=20; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

%% Creating model geometry and mesh
% A box is created with tri-linear hexahedral (hex8) elements using the
% |hexMeshBox| function. The function offers the boundary faces with
% seperate labels for the top, bottom, left, right, front, and back sides.
% As such these can be used to define boundary conditions on the exterior. 

% Create a box with hexahedral elements
cubeDimensions=[sampleWidth sampleThickness sampleHeight]; %Dimensions
cubeElementNumbers=[numElementsWidth numElementsThickness numElementsHeight]; %Number of elements
outputStructType=2; %A structure compatible with mesh view
[meshStruct]=hexMeshBox(cubeDimensions,cubeElementNumbers,outputStructType);

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
FX=Fb(Cb==1,:); %Face orthogonal to X
FY=Fb(Cb==3,:); %Face orthogonal to Y
FZ=Fb(Cb==5,:); %Face orthogonal to Z

indAll=1:1:size(V,1);
bcSupportList_X=indAll(ismember(indAll,FX) & ismember(indAll,FZ)); %Node set part of selected face
bcSupportList_Y=indAll(ismember(indAll,FY) & ismember(indAll,FZ)); %Node set part of selected face
bcSupportList_Z=unique(FZ(:)); %Node set part of selected face

%Prescribed displacement nodes
F_fix=Fb(Cb==6,:); %The current face set
bcPrescribeList=unique(F_fix(:)); %Node set part of selected face

%% 
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists. 

hf=cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,V,'kw','k',0.5);

hl(1)=plotV(V(bcSupportList_X,:),'r.','MarkerSize',markerSize*2);
hl(2)=plotV(V(bcSupportList_Y,:),'g.','MarkerSize',markerSize*2);
hl(3)=plotV(V(bcSupportList_Z,:),'b.','MarkerSize',markerSize);
hl(4)=plotV(V(bcPrescribeList,:),'k.','MarkerSize',markerSize);

legend(hl,{'BC x support','BC y support','BC z support','BC z prescribe'});

axisGeom(gca,fontSize);
camlight headlight; 
drawnow; 

%% DEFINE FIBRE DIRECTIONS

[R,~]=euler2DCM([0,alphaFib,0]);
v_fib=(R*[0 0 1]')';
V_fib=v_fib(ones(size(E,1),1),:);
[a,d]=vectorOrthogonalPair(V_fib);

[VE]=patchCentre(E,V);

%% 
% Visualizing fiber directions

hf=cFigure;
title('Fiber directions','FontSize',fontSize);
hold on;

gpatch(Fb,V,'kw','none',0.25);
[h]=quiverVec(VE,V_fib,1,'k');
[h]=quiverVec(VE,a,0.5,'r');
[h]=quiverVec(VE,d,0.5,'g');

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

%Material section
febio_spec.Material.material{1}.ATTR.type='solid mixture';
febio_spec.Material.material{1}.ATTR.id=1;

febio_spec.Material.material{1}.mat_axis.ATTR.type='user';

febio_spec.Material.material{1}.solid{1}.ATTR.type='Ogden unconstrained';
febio_spec.Material.material{1}.solid{1}.c1=c1;
febio_spec.Material.material{1}.solid{1}.m1=m1;
febio_spec.Material.material{1}.solid{1}.c2=c1;
febio_spec.Material.material{1}.solid{1}.m2=-m1;
febio_spec.Material.material{1}.solid{1}.cp=k;

switch fiberType
    case 1
        febio_spec.Material.material{1}.solid{2}.ATTR.type='fiber-exp-pow';
        febio_spec.Material.material{1}.solid{2}.ksi=ksi;
        febio_spec.Material.material{1}.solid{2}.alpha=alphaPar;
        febio_spec.Material.material{1}.solid{2}.beta=beta;
        febio_spec.Material.material{1}.solid{2}.theta=0;
        febio_spec.Material.material{1}.solid{2}.phi=0;        
    case 2
        febio_spec.Material.material{1}.solid{2}.ATTR.type='ellipsoidal fiber distribution';
        febio_spec.Material.material{1}.solid{2}.ksi=ksi;
        febio_spec.Material.material{1}.solid{2}.beta=beta;
end

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='hex8'; %Element type of this set
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

% -> ElementSets
febio_spec.Geometry.ElementSet{1}.ATTR.name='elementSetTransiso';
febio_spec.Geometry.ElementSet{1}.elem.ATTR.id=(1:size(E,1))';

%MeshData section
% -> ElementData
febio_spec.MeshData.ElementData{1}.ATTR.elem_set=febio_spec.Geometry.ElementSet{1}.ATTR.name;
febio_spec.MeshData.ElementData{1}.ATTR.var='mat_axis';

for q=1:1:size(E,1)
    febio_spec.MeshData.ElementData{1}.elem{q}.ATTR.lid=q;
    febio_spec.MeshData.ElementData{1}.elem{q}.a=a(q,:);
    febio_spec.MeshData.ElementData{1}.elem{q}.d=d(q,:);
end

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

%Output section 
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='s1';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.element_data{1}.VAL=1:size(E,1);

%% Quick viewing of the FEBio input file structure
% The |febView| function can be used to view the xml structure in a MATLAB
% figure window. 

%%
% |febView(febio_spec); %Viewing the febio file|
% febView(febio_spec);

%% Exporting the FEBio input file
% Exporting the febio_spec structure to an FEBio input file is done using
% the |febioStruct2xml| function. 

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
% febView(febioFebFileName); 

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
            
    %%
    % Importing element stress from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),1,1);
    
    %Access data
    E_stress_mat=dataStruct.data;
    
    E_stress_mat(isnan(E_stress_mat))=0;
    
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    [CV]=faceToVertexMeasure(E,V,E_stress_mat(:,:,end));
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('$\sigma_{1}$ [MPa]','Interpreter','Latex')
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
