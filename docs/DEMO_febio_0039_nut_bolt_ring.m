%% DEMO_febio_0039_nut_bolt_ring
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
lineWidth=3;

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');
stlPath=fullfile(defaultFolder,'data','STL');

% STL files for parts
fileName_1=fullfile(stlPath,'M3_nut.stl');
fileName_2=fullfile(stlPath,'M3_bolt.stl');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement

%Define prescribed rotation
numRotations=1;
prescribedRotation_Z=-(pi/2)*numRotations;
prescribedDisplacement_Z=-numRotations*0.5; 

%Material parameter set
c1=1e-3; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
k_factor=1e2; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus
d=1e-9; %Density

% FEA control settings
numTimeSteps=25; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=16; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
symmetric_stiffness=0;
min_residual=1e-20;

%Contact parameters
contactPenalty=30;
laugon=0;
minaug=1;
maxaug=10;
fric_coeff=1; 

%% Import STL file as patch data
[stlStruct] = import_STL(fileName_1);
F_nut=stlStruct.solidFaces{1}; %Faces
V_nut=stlStruct.solidVertices{1}; %Vertices
[F_nut,V_nut]=mergeVertices(F_nut,V_nut); % Merging nodes
V_nut=V_nut*10;
meanV1=mean(V_nut,1);
V_nut=V_nut-meanV1;
% V_nut(:,[1 2])=V_nut(:,[1 2]).*1.01;
V_nut=V_nut+meanV1;

[stlStruct] = import_STL(fileName_2);
F_bolt=stlStruct.solidFaces{1}; %Faces
V_bolt=stlStruct.solidVertices{1}; %Vertices
[F_bolt,V_bolt]=mergeVertices(F_bolt,V_bolt); % Merging nodes
V_bolt=V_bolt*10;

V_nut(:,3)=V_nut(:,3)+4.31-0.25;

%%


pointSpacing=0.3; 

r1=3;
n=round((2*pi*r1)/pointSpacing);
t=linspace(0,2*pi,n);t=t(1:end-1); 
x=r1*cos(t);
y=r1*sin(t);
V1c=[x(:) y(:)];

r2=3.2/2;
n=round((2*pi*r2)/pointSpacing);
t=linspace(0,2*pi,n);t=t(1:end-1); 
x=r2*cos(t);
y=r2*sin(t);
V2c=[x(:) y(:)];

[F1,V1]=regionTriMesh2D({V1c,V2c},pointSpacing,0,0);
F1=fliplr(F1);
V1(:,3)=0;
V1c(:,3)=0;
V2c(:,3)=0;

F2=F1;
F2=fliplr(F2);
V2=V1;
V2(:,3)=V2(:,3)+2; 

V1c2=V1c;
V1c2(:,3)=V1c2(:,3)+2; 

V2c2=V2c;
V2c2(:,3)=V2c2(:,3)+2; 

cPar.closeLoopOpt=1; 
cPar.patchType='tri';
[F3,V3]=polyLoftLinear(V1c,V1c2,cPar);
[F4,V4]=polyLoftLinear(V2c,V2c2,cPar);
F4=fliplr(F4);

[Fw,Vw,Cw]=joinElementSets({F1,F2,F3,F4},{V1,V2,V3,V4});
[Fw,Vw]=mergeVertices(Fw,Vw);

%% Mesh ring
[V_regions]=getInnerPoint(Fw,Vw); % Define region points
V_holes=[]; % Define hole points
[regionA]=tetVolMeanEst(Fw,Vw); %Volume for regular tets

stringOpt='-pq1.2AaY';

inputStruct.stringOpt=stringOpt;
inputStruct.Faces=Fw;
inputStruct.Nodes=Vw;
inputStruct.holePoints=V_holes;
inputStruct.faceBoundaryMarker=Cw; %Face boundary markers
inputStruct.regionPoints=V_regions; %region points
inputStruct.regionA=regionA*10;
inputStruct.minRegionMarker=2; %Minimum region marker

% Mesh model using tetrahedral elements using tetGen 
[meshOutput]=runTetGen(inputStruct); %Run tetGen 

% Access model element and patch data
Fwb=fliplr(meshOutput.facesBoundary);
Vw=meshOutput.nodes;
Cwb=meshOutput.boundaryMarker;
CE=meshOutput.elementMaterialID;
E=meshOutput.elements;

Vw(:,3)=Vw(:,3)+2; 

%% Visualizing mesh using |meshView|, see also |anim8|

meshView(meshOutput);

%% Join geometries

V=[Vw; V_nut; V_bolt];
F_nut=F_nut+size(Vw,1);
F_bolt=F_bolt+size(Vw,1)+size(V_nut,1);

%%
R=euler2DCM([pi 0 0]);
V=V*R; 

%%
% Visualize
cFigure; hold on; 
gpatch(F_nut,V,'rw','k',1);
gpatch(F_bolt,V,'gw','k',1);
gpatch(Fwb,V,Cwb,'k',1);
patchNormPlot(Fwb,V)
colormap(gjet(4)); icolorbar;
axisGeom;
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
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;
febio_spec.Control.max_refs=max_refs;
febio_spec.Control.max_ups=max_ups;
febio_spec.Control.symmetric_stiffness=symmetric_stiffness; 
febio_spec.Control.min_residual=min_residual;

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
febio_spec.Material.material{2}.density=1e-9;
febio_spec.Material.material{2}.center_of_mass=mean(V_nut,1);

febio_spec.Material.material{3}.ATTR.type='rigid body';
febio_spec.Material.material{3}.ATTR.id=3;
febio_spec.Material.material{3}.density=1e-9;
febio_spec.Material.material{3}.center_of_mass=mean(V_bolt,1);

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='tet4'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name='Washer'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=E;

febio_spec.Geometry.Elements{2}.ATTR.type='tri3'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat=2; %material index for this set 
febio_spec.Geometry.Elements{2}.ATTR.name='Nut'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id=((size(E,1)+1):1:(size(E,1)+size(F_nut,1)))'; %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL=F_nut;

febio_spec.Geometry.Elements{3}.ATTR.type='tri3'; %Element type of this set
febio_spec.Geometry.Elements{3}.ATTR.mat=3; %material index for this set 
febio_spec.Geometry.Elements{3}.ATTR.name='Bolt'; %Name of the element set
febio_spec.Geometry.Elements{3}.elem.ATTR.id=((size(E,1)+size(F_nut,1)+1):1:(size(E,1)+size(F_nut,1)+size(F_bolt,1)))'; %Element id's
febio_spec.Geometry.Elements{3}.elem.VAL=F_bolt;

% % -> NodeSets
% febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList';
% febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcSupportList(:);

% -> Surfaces
febio_spec.Geometry.Surface{1}.ATTR.name='contact1_master';
febio_spec.Geometry.Surface{1}.tri3.ATTR.lid=(1:1:size(F_bolt,1))';
febio_spec.Geometry.Surface{1}.tri3.VAL=F_bolt;

febio_spec.Geometry.Surface{2}.ATTR.name='contact_slave';
febio_spec.Geometry.Surface{2}.tri3.ATTR.lid=(1:1:size(Fwb,1))';
febio_spec.Geometry.Surface{2}.tri3.VAL=Fwb;

febio_spec.Geometry.Surface{3}.ATTR.name='contact2_master';
febio_spec.Geometry.Surface{3}.tri3.ATTR.lid=(1:1:size(F_nut,1))';
febio_spec.Geometry.Surface{3}.tri3.VAL=F_nut;

% -> Surface pairs
febio_spec.Geometry.SurfacePair{1}.ATTR.name='Contact1';
febio_spec.Geometry.SurfacePair{1}.master.ATTR.surface=febio_spec.Geometry.Surface{1}.ATTR.name;
febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;

febio_spec.Geometry.SurfacePair{2}.ATTR.name='Contact2';
febio_spec.Geometry.SurfacePair{2}.master.ATTR.surface=febio_spec.Geometry.Surface{3}.ATTR.name;
febio_spec.Geometry.SurfacePair{2}.slave.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;

% %Boundary condition section 
% % -> Fix boundary conditions
% febio_spec.Boundary.fix{1}.ATTR.bc='x';
% febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
% febio_spec.Boundary.fix{2}.ATTR.bc='y';
% febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
% febio_spec.Boundary.fix{3}.ATTR.bc='z';
% febio_spec.Boundary.fix{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;

% -> Prescribed boundary conditions on the rigid body
febio_spec.Boundary.rigid_body{1}.ATTR.mat=2;
febio_spec.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='x';
febio_spec.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='z';
febio_spec.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Rx';
febio_spec.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Ry';
febio_spec.Boundary.rigid_body{1}.fixed{6}.ATTR.bc='Rz';

febio_spec.Boundary.rigid_body{2}.ATTR.mat=3;
febio_spec.Boundary.rigid_body{2}.fixed{1}.ATTR.bc='x';
febio_spec.Boundary.rigid_body{2}.fixed{2}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{2}.fixed{3}.ATTR.bc='z';
febio_spec.Boundary.rigid_body{2}.fixed{4}.ATTR.bc='Rx';
febio_spec.Boundary.rigid_body{2}.fixed{5}.ATTR.bc='Ry';

febio_spec.Boundary.rigid_body{2}.prescribed{1}.ATTR.bc='Rz';
febio_spec.Boundary.rigid_body{2}.prescribed{1}.ATTR.lc=1;
febio_spec.Boundary.rigid_body{2}.prescribed{1}.VAL=prescribedRotation_Z;
febio_spec.Boundary.rigid_body{2}.prescribed{2}.ATTR.bc='z';
febio_spec.Boundary.rigid_body{2}.prescribed{2}.ATTR.lc=1;
febio_spec.Boundary.rigid_body{2}.prescribed{2}.VAL=prescribedDisplacement_Z;

%Contact section
febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
febio_spec.Contact.contact{1}.ATTR.type='sliding-elastic';
febio_spec.Contact.contact{1}.two_pass=1;
febio_spec.Contact.contact{1}.laugon=laugon;
febio_spec.Contact.contact{1}.tolerance=0.2;
febio_spec.Contact.contact{1}.gaptol=0;
febio_spec.Contact.contact{1}.minaug=minaug;
febio_spec.Contact.contact{1}.maxaug=maxaug;
febio_spec.Contact.contact{1}.search_tol=0.01;
febio_spec.Contact.contact{1}.search_radius=0.1;
febio_spec.Contact.contact{1}.symmetric_stiffness=0;
febio_spec.Contact.contact{1}.auto_penalty=1;
febio_spec.Contact.contact{1}.penalty=contactPenalty;
febio_spec.Contact.contact{1}.fric_coeff=fric_coeff;

febio_spec.Contact.contact{2}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{2}.ATTR.name;
febio_spec.Contact.contact{2}.ATTR.type='sliding-elastic';
febio_spec.Contact.contact{2}.two_pass=1;
febio_spec.Contact.contact{2}.laugon=laugon;
febio_spec.Contact.contact{2}.tolerance=0.2;
febio_spec.Contact.contact{2}.gaptol=0;
febio_spec.Contact.contact{2}.minaug=minaug;
febio_spec.Contact.contact{2}.maxaug=maxaug;
febio_spec.Contact.contact{2}.search_tol=0.01;
febio_spec.Contact.contact{2}.search_radius=0.1;
febio_spec.Contact.contact{2}.symmetric_stiffness=0;
febio_spec.Contact.contact{2}.auto_penalty=1;
febio_spec.Contact.contact{2}.penalty=contactPenalty;
febio_spec.Contact.contact{2}.fric_coeff=fric_coeff;


%LoadData section
febio_spec.LoadData.loadcurve{1}.ATTR.id=1;
febio_spec.LoadData.loadcurve{1}.ATTR.type='linear';
febio_spec.LoadData.loadcurve{1}.point.VAL=[0 0; 1 1];

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
febioAnalysis.maxLogCheckTime=3; %Max log file checking time

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
    DN=N_disp_mat(:,:,end);
    DN_magnitude=sqrt(sum(DN.^2,2));
    V_def=V+DN;
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
    X_DEF=V_DEF(:,1,:);
    Y_DEF=V_DEF(:,2,:);
    Z_DEF=V_DEF(:,3,:);
    
    %%
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    % Create basic view and store graphics handle to initiate animation
    [t,p,R]=cart2sph(V_def(:,1),V_def(:,2),V_def(:,3));
    
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1=gpatch(F_bolt,V_def,'kw','none',1); %Add graphics object to animate
    gpatch(F_nut,V_def,'kw','none',0.5); %Add graphics object to animate
    hp2=gpatch(Fwb,V_def,DN_magnitude,'k',1); %Add graphics object to animate
    hp2.FaceColor='interp';
    
    axisGeom(gca,fontSize); 
    colormap(jet(250)); colorbar;caxis([0 max(DN_magnitude)]);
    caxis([0 2.5]);%max(DN_magnitude)]);
    axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);
    camlight headlight;
        
    % Set up animation features
    animStruct.Time=time_mat; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        DN=N_disp_mat(:,:,qt); %Current displacement
        DN_magnitude=sqrt(sum(DN.^2,2)); %Current displacement magnitude
        V_def=V+DN; %Current nodal coordinates
%         [CF]=vertexToFaceMeasure(F,DN_magnitude); %Current color data to use
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1 hp2 hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,V_def,DN_magnitude}; %Property values for to set in order to animate
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
% Copyright (C) 2018  Kevin Mattheus Moerman
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
