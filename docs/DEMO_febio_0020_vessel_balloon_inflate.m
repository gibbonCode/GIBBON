%% DEMO_febio_0020_vessel_balloon_inflate
% Below is a demonstration for:
%
% * Building geometry for a cylindrical vessel with hexahedral elements
% * Defining the boundary conditions
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement results

%% Keywords
%
% * febio_spec version 4.0
% * febio, FEBio
% * vessel, cylinder, balloon
% * prescribed displacement
% * contact, sliding
% * hexahedral elements, hex8
% * tube, cylindrical
% * static, solid, multistep
% * hyperelastic, Ogden
% * displacement logfile
% * stress logfile
% * Nodal displacement mapping
% * Shell thickness specification

%%

clear; close all; clc;

%% Plot settings
fontSize=20;
faceAlpha1=0.8;
markerSize=40;
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
febioLogFileName_strain=[febioFebFileNamePart,'_strain_out.txt']; %Log file name for exporting strain

%Specifying geometry parameters vessel (mm)
pointSpacing=1;
radiusOuter1=6.5/2;
radiusInner1=5.4/2;
radiusOuter2=5.8/2;
radiusInner2=4.7/2;
vesselLength=40;
contactInitialOffset=0.1;

radiusBalloon=min([radiusInner1 radiusInner2])-contactInitialOffset;
pointSpacingBalloon=pointSpacing/2;
balloonExtraLength=pointSpacing/2;

%Define applied bc's
radialDisplacement=1; %radial displacement after touch
radialDisplacementTotal=radialDisplacement+contactInitialOffset; %Total radial displacement

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
contactPenalty=0.5;
laugon=0;
minaug=1;
maxaug=10;
fric_coeff=0;

%Visualization parameters
colorLimits_volumeRatio=[0.99 1.01];

%% Creating model boundary polygons
%

nRad=round((2*pi*max([radiusOuter1 radiusOuter2]))/pointSpacing); %Number of radial steps

t=linspace(0,2*pi,nRad)'; %Angles
t=t(1:end-1); %take away last which equals start
v1_Outer=[-(vesselLength/2)*ones(size(t)) radiusOuter1*sin(t) radiusOuter1*cos(t)]; %Circular coordinates

t=linspace(0,2*pi,nRad)'; %Angles
t=t(1:end-1); %take away last which equals start
v2_Outer=[(vesselLength/2)*ones(size(t)) radiusOuter2*sin(t) radiusOuter2*cos(t)]; %Circular coordinates

t=linspace(0,2*pi,nRad)'; %Angles
t=t(1:end-1); %take away last which equals start
v1_Inner=[-(vesselLength/2)*ones(size(t)) radiusInner1*sin(t) radiusInner1*cos(t)]; %Circular coordinates

t=linspace(0,2*pi,nRad)'; %Angles
t=t(1:end-1); %take away last which equals start
v2_Inner=[(vesselLength/2)*ones(size(t)) radiusInner2*sin(t) radiusInner2*cos(t)]; %Circular coordinates

t=linspace(0,2*pi,nRad)'; %Angles
t=t(1:end-1); %take away last which equals start
v1_balloon=[-((vesselLength/2)+balloonExtraLength)*ones(size(t)) radiusBalloon*sin(t) radiusBalloon*cos(t)]; %Circular coordinates
v2_balloon=[ ((vesselLength/2)+balloonExtraLength)*ones(size(t)) radiusBalloon*sin(t) radiusBalloon*cos(t)]; %Circular coordinates

%%
% Plotting model boundary polygons

cFigure;
hold on;
title('Model boundary polygons','FontSize',fontSize);
plotV(v1_Outer,'r.-')
plotV(v1_Inner,'g.-')
plotV(v2_Outer,'b.-')
plotV(v2_Inner,'y.-')
plotV(v1_balloon,'c.-')
plotV(v2_balloon,'k.-')
axisGeom(gca,fontSize);

drawnow;

%% Creating model boundary surfaces

controlStructLoft.numSteps=ceil(vesselLength./pointSpacing);
controlStructLoft.closeLoopOpt=1;
controlStructLoft.patchType='quad';

%Meshing outer surface
[F1,V1,indStart1]=polyLoftLinear(v1_Outer,v2_Outer,controlStructLoft);

%Meshing inner surface
[F2,V2,indStart2]=polyLoftLinear(v1_Inner,v2_Inner,controlStructLoft);

%Compose hexahedral elements
indStart2=indStart2+size(V1,1);
F2=F2+size(V1,1);
E=[F1 F2];
V=[V1;V2];
[FE]=element2patch(E,[],'hex8');
indBoundary=tesBoundary(FE);
Fb=FE(indBoundary,:);

%Meshing balloon surface
[Fs,Vs]=polyLoftLinear(v1_balloon,v2_balloon,controlStructLoft);
Fs=fliplr(Fs); %Invert orientation
[Fs,Vs]=subQuad(Fs,Vs,1);

%%
% Plotting model boundary surfaces

cFigure;
hold on;
title('Model boundary surfaces','FontSize',fontSize);

gpatch(FE,V,'rw','k',0.5);
patchNormPlot(FE,V);

gpatch(F1,V,'bw');
patchNormPlot(F1,V);

gpatch(F2,V,'gw');
patchNormPlot(F2,V);

gpatch(Fs,Vs,'kw');
patchNormPlot(Fs,Vs);

colormap(gca,gjet(4));
icolorbar;

axisGeom(gca,fontSize);
camlight headlight;
drawnow;

%% Joining node sets

Fs=Fs+size(V,1);
V=[V;Vs];

%% Defining the boundary conditions
% The visualization of the model boundary shows colors for each side of the
% cube. These labels can be used to define boundary conditions.

%Define X supported node set
bcSupportList_X=unique([indStart1(:);indStart2(:)]); %Node set part of selected face

% %Define Y supported node set
% bcSupportList_Y=unique(Fb(Cb==4,:)); %Node set part of selected face

%Radial expansion prescribed displacement
bcPrescribeList=(size(V,1)-size(Vs,1)+1):size(V,1);

%Define presribed displacements
[t,r,z]  = cart2pol(V(bcPrescribeList,2),V(bcPrescribeList,3),V(bcPrescribeList,1));
[y2,z2,x2] = pol2cart(t,r+radialDisplacementTotal,z);
Ux=x2-V(bcPrescribeList,1);
Uy=y2-V(bcPrescribeList,2);
Uz=z2-V(bcPrescribeList,3);

%%
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists.

cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(FE,V,'kw','k',0.5);
gpatch(Fs,V,'kw','k',0.5);

hl(1)=plotV(V(bcSupportList_X,:),'k.','MarkerSize',markerSize);
hl(2)=gpatch(Fs,V,'rw','r',0.5);

legend(hl,{'BC X support','BC prescribe'});

axisGeom(gca,fontSize);
camlight headlight;
drawnow;

%% Define contact

F_contact_primary=Fs;
F_contact_secondary=F2;

%%
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists.

hf=cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(FE,V,'kw','none',0.5);

h1(1)=gpatch(F_contact_primary,V,'gw','g',0.5);
patchNormPlot(F_contact_primary,V);
h1(2)=gpatch(F_contact_secondary,V,'rw','r',0.5);
patchNormPlot(F_contact_secondary,V);

legend(hl,{'Primary surface','Secondary surface'});

axisGeom(gca,fontSize);
camlight headlight;
drawnow;

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
febio_spec.Material.material{1}.c2=c1;
febio_spec.Material.material{1}.m2=-m1;
febio_spec.Material.material{1}.k=k;

materialName2='Material2';
febio_spec.Material.material{2}.ATTR.name=materialName2;
febio_spec.Material.material{2}.ATTR.type='Ogden';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.c1=c1;
febio_spec.Material.material{2}.m1=m1;
febio_spec.Material.material{2}.c2=c1;
febio_spec.Material.material{2}.m2=-m1;
febio_spec.Material.material{2}.k=k;

%Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type 
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix

partName2='Part2';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='quad4'; %Element type 
febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E,1)+(1:1:size(Fs,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=Fs; %The element matrix

% -> Surfaces
surfaceName1='contactSurface1';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(F_contact_primary,1))';
febio_spec.Mesh.Surface{1}.quad4.VAL=F_contact_primary;

surfaceName2='contactSurface2';
febio_spec.Mesh.Surface{2}.ATTR.name=surfaceName2;
febio_spec.Mesh.Surface{2}.quad4.ATTR.id=(1:1:size(F_contact_secondary,1))';
febio_spec.Mesh.Surface{2}.quad4.VAL=F_contact_secondary;

% -> Surface pairs
febio_spec.Mesh.SurfacePair{1}.ATTR.name='Contact1';
febio_spec.Mesh.SurfacePair{1}.primary=surfaceName1;
febio_spec.Mesh.SurfacePair{1}.secondary=surfaceName2;

% -> NodeSets
nodeSetName1='bcSupportList_X';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcSupportList_X);

nodeSetName2='bcPrescribeList';
febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.VAL=mrow(bcPrescribeList);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

febio_spec.MeshDomains.ShellDomain.ATTR.name=partName2;
febio_spec.MeshDomains.ShellDomain.ATTR.mat=materialName2;
febio_spec.MeshDomains.ShellDomain.shell_thickness=0.01;

%-> Node data
meshDataName_Ux='DisplacementData_X';           
febio_spec.MeshData.NodeData{1}.ATTR.name=meshDataName_Ux;
febio_spec.MeshData.NodeData{1}.ATTR.node_set=nodeSetName2;
febio_spec.MeshData.NodeData{1}.ATTR.data_type='scalar';
febio_spec.MeshData.NodeData{1}.node.ATTR.lid=(1:numel(bcPrescribeList))';
febio_spec.MeshData.NodeData{1}.node.VAL=Ux(:);

meshDataName_Uy='DisplacementData_Y';           
febio_spec.MeshData.NodeData{2}.ATTR.name=meshDataName_Uy;
febio_spec.MeshData.NodeData{2}.ATTR.node_set=nodeSetName2;
febio_spec.MeshData.NodeData{2}.ATTR.data_type='scalar';
febio_spec.MeshData.NodeData{2}.node.ATTR.lid=(1:numel(bcPrescribeList))';
febio_spec.MeshData.NodeData{2}.node.VAL=Uy(:);

meshDataName_Uz='DisplacementData_Z';           
febio_spec.MeshData.NodeData{3}.ATTR.name=meshDataName_Uz;
febio_spec.MeshData.NodeData{3}.ATTR.node_set=nodeSetName2;
febio_spec.MeshData.NodeData{3}.ATTR.data_type='scalar';
febio_spec.MeshData.NodeData{3}.node.ATTR.lid=(1:numel(bcPrescribeList))';
febio_spec.MeshData.NodeData{3}.node.VAL=Uz(:);

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_x';
febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.x_dof=1;
febio_spec.Boundary.bc{1}.y_dof=0;
febio_spec.Boundary.bc{1}.z_dof=0;

% -> Prescribe boundary conditions
febio_spec.Boundary.bc{2}.ATTR.name='prescibed_displacement_x';
febio_spec.Boundary.bc{2}.ATTR.type='prescribed displacement';
febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{2}.dof='x';
febio_spec.Boundary.bc{2}.value.ATTR.lc=1;
febio_spec.Boundary.bc{2}.value.ATTR.type='map';
febio_spec.Boundary.bc{2}.value.VAL=meshDataName_Ux;
febio_spec.Boundary.bc{2}.relative=0;

febio_spec.Boundary.bc{3}.ATTR.name='prescibed_displacement_y';
febio_spec.Boundary.bc{3}.ATTR.type='prescribed displacement';
febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{3}.dof='y';
febio_spec.Boundary.bc{3}.value.ATTR.lc=1;
febio_spec.Boundary.bc{3}.value.ATTR.type='map';
febio_spec.Boundary.bc{3}.value.VAL=meshDataName_Uy;
febio_spec.Boundary.bc{3}.relative=0;

febio_spec.Boundary.bc{4}.ATTR.name='prescibed_displacement_z';
febio_spec.Boundary.bc{4}.ATTR.type='prescribed displacement';
febio_spec.Boundary.bc{4}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{4}.dof='z';
febio_spec.Boundary.bc{4}.value.ATTR.lc=1;
febio_spec.Boundary.bc{4}.value.ATTR.type='map';
febio_spec.Boundary.bc{4}.value.VAL=meshDataName_Uz;
febio_spec.Boundary.bc{4}.relative=0;

%Contact section
febio_spec.Contact.contact{1}.ATTR.type='sliding-elastic';
febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Mesh.SurfacePair{1}.ATTR.name;
febio_spec.Contact.contact{1}.two_pass=1;
febio_spec.Contact.contact{1}.laugon=laugon;
febio_spec.Contact.contact{1}.tolerance=0.2;
febio_spec.Contact.contact{1}.gaptol=0;
febio_spec.Contact.contact{1}.minaug=minaug;
febio_spec.Contact.contact{1}.maxaug=maxaug;
febio_spec.Contact.contact{1}.search_tol=0.01;
febio_spec.Contact.contact{1}.search_radius=0.1*sqrt(sum((max(V,[],1)-min(V,[],1)).^2,2));
febio_spec.Contact.contact{1}.symmetric_stiffness=0;
febio_spec.Contact.contact{1}.auto_penalty=1;
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

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_strain;
febio_spec.Output.logfile.element_data{1}.ATTR.data='E1';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

febio_spec.Output.plotfile.compression=0;

%% Quick viewing of the FEBio input file structure
% The |febView| function can be used to view the xml structure in a MATLAB
% figure window.

%%
% |febView(febio_spec); %Viewing the febio file|
% febView(febio_spec)

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
febioAnalysis.runMode=runMode;%'internal';

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
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    DN_magnitude=sqrt(sum(N_disp_mat(:,:,end).^2,2)); %Current displacement magnitude
        
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('Displacement magnitude [mm]','Interpreter','Latex')
    hp=gpatch(Fb,V_DEF(:,:,end),DN_magnitude,'k',1); %Add graphics object to animate
    hp.FaceColor='interp';
    hp2=gpatch(Fs,V,'w','k',1); %A static graphics object
    
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
        animStruct.Handles{qt}=[hp hp hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),DN_magnitude,V_DEF(:,:,qt)}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;
            
    %%
    % Importing element stress from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_strain),0,1);
    
    %Access data
    E_strain_mat=dataStruct.data;
    E_strain_mat(isnan(E_strain_mat))=0;
    
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    [CV]=faceToVertexMeasure(E,V,E_strain_mat(:,:,end));
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('$\varepsilon_{1}$ [.]','Interpreter','Latex')
    hp=gpatch(Fb,V_DEF(:,:,end),CV,'k',1); %Add graphics object to animate
    hp.FaceColor='interp';
    hp2=gpatch(Fs,V,'w','k',1); %A static graphics object
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([min(E_strain_mat(:)) max(E_strain_mat(:))]);    
    axis(axisLim(V_DEF)); %Set axis limits statically    
    camlight headlight;        
        
    % Set up animation features
    animStruct.Time=timeVec; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        
        [CV]=faceToVertexMeasure(E,V,E_strain_mat(:,:,qt));
        
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
