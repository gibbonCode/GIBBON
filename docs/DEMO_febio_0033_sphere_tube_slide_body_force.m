%% DEMO_febio_0033_sphere_tube_slide_body_force
% Below is a demonstration for:
% 
% * Building geometry for a spherical blob with tetrahedral elements
% which is being aspirated into a tube. 
% This demo consists off:
% * Defining the boundary conditions 
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement results

%% Keywords
%
% * febio_spec version 4.0
% * febio, FEBio
% * indentation
% * contact, sliding, friction
% * rigid body constraints
% * hexahedral elements, hex8
% * quadrilateral elements, quad4
% * shell elements
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

% Sphere parameters
sphereRadius=5; 
pointSpacing=0.4; 

% Ground plate parameters
tubeRadius=0.8*sphereRadius; 
inletRadius=tubeRadius/3;
tubeLength=sphereRadius*6; 

% Material parameter set
materialType=1; 
c1=0.1e-3; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
k_factor=10; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus
g1=2; %Viscoelastic QLV proportional coefficient
t1=150; %Viscoelastic QLV time coefficient

%Setting material density
materialDensity=1e-9; %Density

% FEA control settings
timeTotal=1; %Total simulation time
numTimeSteps=40; %Number of time steps desired
max_refs=20; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(timeTotal/numTimeSteps)/100; %Minimum time step size
dtmax=timeTotal/(numTimeSteps); %Maximum time step size
symmetric_stiffness=0;
min_residual=1e-20;
analysisType='DYNAMIC';

runMode='external';% 'internal' or 'external'

%Contact parameters
contactPenalty=50;
fric_coeff=0.25; 
laugon=0;
minaug=1;
maxaug=10;

%Specifying load
sphereVolume=4/3*(pi*sphereRadius^3); %Sphere Volume in mm^3
sphereMass=sphereVolume.*materialDensity; %Sphere mass in tonne
sphereSectionArea=pi*sphereRadius^2;
bodyLoadMagnitude=(9.81*sphereMass*1000)*6500000;

forceBodyLoad=sphereMass.*bodyLoadMagnitude;
stressBodyLoad=forceBodyLoad/sphereSectionArea;

%% Creating model geometry and mesh
% Creating a solid hexahedral mesh sphere

%Control settings
optionStruct.sphereRadius=sphereRadius;
optionStruct.coreRadius=sphereRadius.*0.75;
optionStruct.numElementsCore=ceil((sphereRadius/2)/pointSpacing); 
optionStruct.numElementsMantel=ceil((sphereRadius-optionStruct.coreRadius)/(2*pointSpacing)); 
optionStruct.makeHollow=0;
optionStruct.outputStructType=2;

%Creating sphere
[meshOutput]=hexMeshSphere(optionStruct);

% Access model element and patch data
Fb_blob=meshOutput.facesBoundary;
Cb_blob=meshOutput.boundaryMarker;
V_blob=meshOutput.nodes;
E_blob=meshOutput.elements;

%%
% Visualize mesh

hFig=cFigure; 
subplot(1,2,1); hold on;
title('Boundary surfaces','FontSize',fontSize);
gpatch(Fb_blob,V_blob,Cb_blob,'k',0.5);
% patchNormPlot(Fb,V);
axisGeom(gca,fontSize);
colormap(gjet); icolorbar;
camlight headlight; 

hs=subplot(1,2,2); hold on; 
title('Cut view of solid mesh','FontSize',fontSize);
optionStruct.hFig=[hFig hs];
gpatch(Fb_blob,V_blob,'kw','none',0.25);
meshView(meshOutput,optionStruct);
axisGeom(gca,fontSize);
drawnow; 

%% Creating tube model
% 
pointSpacingTube=mean(patchEdgeLengths(Fb_blob,V_blob))/2;
t=linspace(-0.1*pi,pi,100);
x=inletRadius*sin(t);
y=inletRadius*cos(t);
V_curve_tube=[x(:) y(:) zeros(size(x(:)))];
V_curve_tube(:,1)=V_curve_tube(:,1)-inletRadius;
V_curve_tube(:,2)=V_curve_tube(:,2)+inletRadius+tubeRadius;
V_curve_tube(end+1,:)=[-tubeLength tubeRadius 0];
V_curve_tube(end+1,:)=[-tubeLength tubeRadius/6 0];
nResample=ceil(max(pathLength(V_curve_tube))./pointSpacingTube);
V_curve_tube=evenlySampleCurve(V_curve_tube,nResample,'pchip',0);

cPar.closeLoopOpt=1;
cPar.numSteps=[]; %If empty the number of steps is derived from point spacing of input curve
cPar.w=[1 0 0];
[F_tube,V_tube]=polyRevolve(V_curve_tube,cPar);
V_tube(:,1)=V_tube(:,1)-sphereRadius;

d=0.1; %Tube placement increment/tolerance
while 1
    R=sqrt(sum(V_tube.^2,2)); %Get current radii    
    if min(R)<sphereRadius %If the tube is in the sphere space
        break
    else %If not in sphere space
        V_tube(:,1)=V_tube(:,1)+d; %Shift a bit
    end
end
center_of_mass_tube=mean(V_tube,1);

%%
% Visualization

cFigure; hold on;
gtitle('The surface meshes',fontSize);
gpatch(Fb_blob,V_blob,'kw','none',0.5);
gpatch(F_tube,V_tube,'bw','none',0.5);
% patchNormPlot(F_tube,V_tube);
axisGeom(gca,fontSize);
camlight headlight; 
drawnow; 

%% Join model node sets

V=[V_blob; V_tube; ];
F_tube=F_tube+size(V_blob,1);

%%
% Visualizing model

cFigure; hold on;
gtitle('Model components',fontSize);
hl(1)=gpatch(Fb_blob,V,'rw','k',0.8);
hl(2)=gpatch(F_tube,V,'bw','k',0.8);
legend(hl,{'Blob','Tube'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 
drawnow; 

%% Get contact surfaces
%

F_contact_secondary=Fb_blob;

%%
% Visualize contact surfaces

cFigure; hold on;
title('The contact pair','fontsize',fontSize);
hl(1)=gpatch(F_tube,V,'rw','k',1);
patchNormPlot(F_tube,V);
hl(2)=gpatch(F_contact_secondary,V,'gw','k',1);
patchNormPlot(F_contact_secondary,V);
legend(hl,{'Secondary','Primary'}); clear hl;
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
febio_spec.Control.analysis=analysisType;
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=timeTotal/numTimeSteps;
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
switch materialType
    case 0
        febio_spec.Material.material{1}.ATTR.type='Ogden unconstrained';
        febio_spec.Material.material{1}.ATTR.id=1;
        febio_spec.Material.material{1}.c1=c1;
        febio_spec.Material.material{1}.m1=m1;
        febio_spec.Material.material{1}.c2=c1;
        febio_spec.Material.material{1}.m2=-m1;
        febio_spec.Material.material{1}.cp=k;
        febio_spec.Material.material{1}.density=materialDensity;
    case 1
        %Viscoelastic part
        febio_spec.Material.material{1}.ATTR.type='viscoelastic';        
        febio_spec.Material.material{1}.ATTR.id=1;
        febio_spec.Material.material{1}.g1=g1;
        febio_spec.Material.material{1}.t1=t1;
        febio_spec.Material.material{1}.density=materialDensity;
        
        %Elastic part
        febio_spec.Material.material{1}.elastic{1}.ATTR.type='Ogden unconstrained';
        febio_spec.Material.material{1}.elastic{1}.c1=c1;
        febio_spec.Material.material{1}.elastic{1}.m1=m1;
        febio_spec.Material.material{1}.elastic{1}.c2=c1;
        febio_spec.Material.material{1}.elastic{1}.m2=-m1;
        febio_spec.Material.material{1}.elastic{1}.cp=k;
        febio_spec.Material.material{1}.elastic{1}.density=materialDensity;
end

materialName2='Material2';
febio_spec.Material.material{2}.ATTR.name=materialName2;
febio_spec.Material.material{2}.ATTR.type='rigid body';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.density=1;
febio_spec.Material.material{2}.center_of_mass=center_of_mass_tube;

%Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type 
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E_blob,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E_blob; %The element matrix

partName2='Part2';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='quad4'; %Element type 
febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E_blob,1)+(1:1:size(F_tube,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=F_tube; %The element matrix

% -> Surfaces
surfaceName1='contactSurface1';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(F_tube,1))';
febio_spec.Mesh.Surface{1}.quad4.VAL=F_tube;

surfaceName2='contactSurface2';
febio_spec.Mesh.Surface{2}.ATTR.name=surfaceName2;
febio_spec.Mesh.Surface{2}.quad4.ATTR.id=(1:1:size(F_contact_secondary,1))';
febio_spec.Mesh.Surface{2}.quad4.VAL=F_contact_secondary;

% -> Surface pairs
contactPairName='Contact1';
febio_spec.Mesh.SurfacePair{1}.ATTR.name=contactPairName;
febio_spec.Mesh.SurfacePair{1}.primary=surfaceName2;
febio_spec.Mesh.SurfacePair{1}.secondary=surfaceName1;

%MeshDomains section
febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

febio_spec.MeshDomains.ShellDomain.ATTR.name=partName2;
febio_spec.MeshDomains.ShellDomain.ATTR.mat=materialName2;

%Loads section
% -> Body load
febio_spec.Loads.body_load.ATTR.type='const';
febio_spec.Loads.body_load.x.ATTR.lc=1;
febio_spec.Loads.body_load.x.VAL=bodyLoadMagnitude;
febio_spec.Loads.body_load.y.ATTR.lc=1;
febio_spec.Loads.body_load.y.VAL=0;
febio_spec.Loads.body_load.z.ATTR.lc=1;
febio_spec.Loads.body_load.z.VAL=0;

%Rigid section 
% ->Rigid body fix boundary conditions
febio_spec.Rigid.rigid_bc{1}.ATTR.name='RigidFix_1';
febio_spec.Rigid.rigid_bc{1}.ATTR.type='rigid_fixed';
febio_spec.Rigid.rigid_bc{1}.rb=2;
febio_spec.Rigid.rigid_bc{1}.Rx_dof=1;
febio_spec.Rigid.rigid_bc{1}.Ry_dof=1;
febio_spec.Rigid.rigid_bc{1}.Rz_dof=1;
febio_spec.Rigid.rigid_bc{1}.Ru_dof=1;
febio_spec.Rigid.rigid_bc{1}.Rv_dof=1;
febio_spec.Rigid.rigid_bc{1}.Rw_dof=1;

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
febio_spec.Contact.contact{1}.penalty=contactPenalty;
febio_spec.Contact.contact{1}.fric_coeff=fric_coeff;

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.name='LC1';
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
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

% Plotfile section
febio_spec.Output.plotfile.compression=0;

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
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    DN_magnitude=sqrt(sum(N_disp_mat(:,:,end).^2,2)); %Current displacement magnitude
        
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; hold on;
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp=gpatch(Fb_blob,V_DEF(:,:,end),DN_magnitude,'k',1); %Add graphics object to animate
    hp.FaceColor='interp';
    gpatch(F_tube,V,'w','none',0.5); %Add graphics object to animate
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]/3);
    axis(axisLim(V_DEF)); %Set axis limits statically  
    camlight headlight;
    drawnow; 
    
    % Set up animation features
    animStruct.Time=timeVec; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments
        DN_magnitude=sqrt(sum(N_disp_mat(:,:,qt).^2,2)); %Current displacement magnitue
        
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
