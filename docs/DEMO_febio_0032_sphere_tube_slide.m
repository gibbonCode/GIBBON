%% DEMO_febio_0032_sphere_tube_slide
% Below is a demonstration for:
% 
% * Building geometry for a spherical blob with tetrahedral elements
% which is being pushed into a tube. 
% This demo consists off:
% * Defining the boundary conditions 
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement results

%% Keywords
%
% * febio_spec version 2.5
% * febio, FEBio
% * indentation
% * contact, sliding, friction
% * rigid body constraints
% * tetrahedral elements, tet4
% * triangular elements, tri3
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
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
% febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

% Sphere parameters
sphereRadius1=1; 
nRefine1=2;
sphereRadius2=sphereRadius1/4; 
nRefine2=1;

% Ground plate parameters
tubeRadius=0.65*sphereRadius1; 
inletRadius=tubeRadius/3;
tubeLength=4; 

% Material parameter set
c1=1e-3; %Shear-modulus-like parameter
m1=8; %Material parameter setting degree of non-linearity
k_factor=10; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus

% FEA control settings
numTimeSteps=25; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=15; %Optimum number of iterations
max_retries=10; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/(numTimeSteps); %Maximum time step size
symmetric_stiffness=0;

contactPenalty=5;

displacement_X=-2.5;

%% Creating model geometry and mesh
% 

% Create hemi-sphere mesh
[F_blob1,V_blob1]=geoSphere(nRefine1,sphereRadius1);
[F_blob2,V_blob2]=geoSphere(nRefine2,sphereRadius2);
[F_blob,V_blob,C_blob]=joinElementSets({F_blob1,F_blob2},{V_blob1,V_blob2});
pointSpacingBlob=mean(patchEdgeLengths(F_blob,V_blob));

%%
% Visualize hemi-sphere surface

cFigure; hold on;
gtitle('The sphere surface mesh',fontSize);
gpatch(F_blob,V_blob,C_blob,'k',0.85);
patchNormPlot(F_blob,V_blob);
axisGeom(gca,fontSize);
colormap(gjet); icolorbar;
camlight headlight; 
drawnow; 

%% 
% Using tetgen to create a tetrahedral mesh

% Tetgen input structure
inputStruct.stringOpt='-pq1.2AaY';
inputStruct.Faces=F_blob;
inputStruct.Nodes=V_blob;
inputStruct.holePoints=[0 0 0];
inputStruct.faceBoundaryMarker=C_blob; %Face boundary markers
inputStruct.regionPoints=getInnerPoint(F_blob,V_blob); %region points
inputStruct.regionA=tetVolMeanEst(F_blob,V_blob);
inputStruct.minRegionMarker=2; %Minimum region marker
 
% Create tetrahedral mesh using tetGen 
[meshOutput]=runTetGen(inputStruct); %Run tetGen 
 
% Access model element and patch data
Fb_blob=fliplr(meshOutput.facesBoundary);
Cb_blob=meshOutput.boundaryMarker;
V_blob=meshOutput.nodes;
E_blob=meshOutput.elements;

%%
% Visualize blob mesh

hFig=cFigure; 
subplot(1,2,1); hold on;
gpatch(Fb_blob,V_blob,Cb_blob,'k',0.8);
patchNormPlot(Fb_blob,V_blob);
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
pointSpacingTube=pointSpacingBlob/2;
t=linspace(-0.1*pi,pi,100);
x=inletRadius*sin(t);
y=inletRadius*cos(t);
V_curve_tube=[x(:) y(:) zeros(size(x(:)))];
V_curve_tube(:,1)=V_curve_tube(:,1)-inletRadius;
V_curve_tube(:,2)=V_curve_tube(:,2)+inletRadius+tubeRadius;
V_curve_tube(end+1,:)=[-tubeLength tubeRadius 0];
nResample=ceil(max(pathLength(V_curve_tube))./pointSpacingTube);
V_curve_tube=evenlySampleCurve(V_curve_tube,nResample,'pchip',0);

cPar.closeLoopOpt=1;
cPar.numSteps=[]; %If empty the number of steps is derived from point spacing of input curve
cPar.w=[1 0 0];
[F_tube,V_tube]=polyRevolve(V_curve_tube,cPar);
V_tube(:,1)=V_tube(:,1)-sphereRadius1;

c=1;
while 1
    [D,indMin]=minDist(V_tube,V_blob);
    [~,indMinMin]=min(D);
    d=V_tube(indMinMin,1)-V_blob(indMin(indMinMin),1);
    V_tube(:,1)=V_tube(:,1)-d(1);
    if c>1
        if abs(d-dp)<0.001
            break
        end
    end
    c=c+1;
    dp=d;
end
center_of_mass_tube=mean(V_tube,1);
%
% Visualizing plate mesh

cFigure; hold on;
gtitle('The plate surface mesh',fontSize);
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

F_contact_blob=Fb_blob(Cb_blob==1,:);

%%
% Visualize contact surfaces

cFigure; hold on;
title('Tube blob contact pair','fontsize',fontSize);
hl(1)=gpatch(F_tube,V,'rw','k',1);
patchNormPlot(F_tube,V);
hl(2)=gpatch(F_contact_blob,V,'gw','k',1);
patchNormPlot(F_contact_blob,V);
legend(hl,{'Master','Slave'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 
drawnow; 

%%
% Visualize contact surfaces

bcPrescribeSet=unique(Fb_blob(Cb_blob==2,:));

cFigure; hold on;
title('Tube blob contact pair','fontsize',fontSize);
hl(1)=gpatch(F_tube,V,'rw','none',0.25);
hl(2)=gpatch(Fb_blob,V,'gw','none',0.25);
hl(3)=plotV(V(bcPrescribeSet,:),'k.','MarkerSize',15);
legend(hl,{'Tube','Blob','Prescribed nodes'}); clear hl;
axisGeom(gca,fontSize);

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
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;
febio_spec.Control.max_refs=max_refs;
febio_spec.Control.max_ups=max_ups;
febio_spec.Control.symmetric_stiffness=symmetric_stiffness; 

%Material section
febio_spec.Material.material{1}.ATTR.type='Ogden unconstrained';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c1=c1;
febio_spec.Material.material{1}.m1=m1;
febio_spec.Material.material{1}.c2=c1;
febio_spec.Material.material{1}.m2=-m1;
febio_spec.Material.material{1}.cp=k;

febio_spec.Material.material{2}.ATTR.type='rigid body';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.density=1;
febio_spec.Material.material{2}.center_of_mass=center_of_mass_tube;

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='tet4'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name='Blob'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(E_blob,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=E_blob;

febio_spec.Geometry.Elements{2}.ATTR.type='quad4'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat=2; %material index for this set 
febio_spec.Geometry.Elements{2}.ATTR.name='Tube'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id=size(E_blob,1)+(1:1:size(F_tube,1))'; %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL=F_tube;

% -> NodeSets
febio_spec.Geometry.NodeSet{1}.ATTR.name='bcMove';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcPrescribeSet(:);

% % -> Surfaces
febio_spec.Geometry.Surface{1}.ATTR.name='contact_master1';
febio_spec.Geometry.Surface{1}.quad4.ATTR.lid=(1:1:size(F_tube,1))';
febio_spec.Geometry.Surface{1}.quad4.VAL=F_tube;

febio_spec.Geometry.Surface{2}.ATTR.name='contact_slave1';
febio_spec.Geometry.Surface{2}.tri3.ATTR.lid=(1:1:size(F_contact_blob,1))';
febio_spec.Geometry.Surface{2}.tri3.VAL=F_contact_blob;

% -> Surface pairs
febio_spec.Geometry.SurfacePair{1}.ATTR.name='Contact1_tube_blob';
febio_spec.Geometry.SurfacePair{1}.master.ATTR.surface=febio_spec.Geometry.Surface{1}.ATTR.name;
febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;

%Boundary condition section 
% -> Fix boundary conditions
% febio_spec.Boundary.fix{1}.ATTR.bc='x';
% febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{1}.ATTR.bc='y';
febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{2}.ATTR.bc='z';
febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
% 
febio_spec.Boundary.prescribe{1}.ATTR.bc='x';
febio_spec.Boundary.prescribe{1}.ATTR.relative=1;
febio_spec.Boundary.prescribe{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.prescribe{1}.scale.ATTR.lc=1;
febio_spec.Boundary.prescribe{1}.scale.VAL=1;
febio_spec.Boundary.prescribe{1}.relative=1;
febio_spec.Boundary.prescribe{1}.value=displacement_X;

% -> Prescribed boundary conditions on the rigid body
febio_spec.Boundary.rigid_body{1}.ATTR.mat=2;
febio_spec.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='x';
febio_spec.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='z';
febio_spec.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Rx';
febio_spec.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Ry';
febio_spec.Boundary.rigid_body{1}.fixed{6}.ATTR.bc='Rz';

%Contact section
% -> Contact 1
febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
febio_spec.Contact.contact{1}.ATTR.type='sliding-elastic';
febio_spec.Contact.contact{1}.two_pass=1;
febio_spec.Contact.contact{1}.laugon=1;
febio_spec.Contact.contact{1}.tolerance=0.2;
febio_spec.Contact.contact{1}.gaptol=0;
febio_spec.Contact.contact{1}.minaug=1;
febio_spec.Contact.contact{1}.maxaug=10;
febio_spec.Contact.contact{1}.search_tol=0.01;
febio_spec.Contact.contact{1}.search_radius=0.1;
febio_spec.Contact.contact{1}.symmetric_stiffness=0;
febio_spec.Contact.contact{1}.auto_penalty=1;
febio_spec.Contact.contact{1}.penalty=contactPenalty;
febio_spec.Contact.contact{1}.fric_coeff=0.1;

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

if 1%runFlag==1 %i.e. a succesful run
    
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
