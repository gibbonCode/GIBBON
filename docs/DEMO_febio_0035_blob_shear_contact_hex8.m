%% DEMO_febio_0035_blob_shear_contact_hex8
% Below is a demonstration for:
% 
% * Building geometry for a hemi-spherical blob with tetrahedral elements
% which is being sheared by a rigid wall. 
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
% * contact, sliding, sticky, friction
% * rigid body constraints
% * tetrahedral elements, tet4
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

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
% febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

% Hemi-sphere parameters
hemiSphereRadius=1;
numElementsMantel=5;
smoothEdge=1; 

% Ground plate parameters
plateRadius=2*hemiSphereRadius; 

% Probe parameters
probeWidth=3*hemiSphereRadius; 
filletProbe=0.5; %Fillet radius

% Define probe displacement
probeDisplacement=hemiSphereRadius*2; 
probeOverlapFactor=0.4;
probeLength=hemiSphereRadius*3; 

% Material parameter set
c1=1e-3; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
k_factor=100; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus

% FEA control settings
numTimeSteps=20;
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=10; %Optimum number of iterations
max_retries=25; %Maximum number of retires
symmetric_stiffness=0;
min_residual=1e-20;
step_size=1/numTimeSteps;
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/(numTimeSteps); %Maximum time step size

%Contact parameters
contactPenalty=5;
laugon=0;
minaug=1;
maxaug=10;
fric_coeff=0.01; 
max_traction=0;  

%% Creating model geometry and mesh
% 

%Control settings
optionStruct.sphereRadius=hemiSphereRadius;
optionStruct.coreRadius=optionStruct.sphereRadius/2;
optionStruct.numElementsMantel=numElementsMantel;
optionStruct.numElementsCore=optionStruct.numElementsMantel*2;
optionStruct.outputStructType=2;
optionStruct.makeHollow=0;
optionStruct.cParSmooth.n=25;

% %Creating sphere
[meshStruct]=hexMeshHemiSphere(optionStruct);

% Access model element and patch data
Fb_blob=meshStruct.facesBoundary;
Cb_blob=meshStruct.boundaryMarker;
V_blob=meshStruct.nodes;
E_blob=meshStruct.elements;

F_blob=element2patch(E_blob);

%%

pointSpacingBlob=mean(patchEdgeLengths(Fb_blob,V_blob));

%Smoothen edges
if smoothEdge==1
    %Get rigid region
    ind=1:1:size(V_blob,1); %Indices for all nodes
    indRigid1=find(ismember(ind,Fb_blob(Cb_blob==2,:)) & ~ismember(ind,Fb_blob(Cb_blob==1,:))); %Indices for new bottom surface nodes
    indRigid2=find(ismember(ind,Fb_blob(Cb_blob==1,:)) & ~ismember(ind,Fb_blob(Cb_blob==2,:))); %Indices for new bottom surface nodes
    indRigid=[indRigid1(:); indRigid2(:);];
    
    %Smoothing
    cPar.Method='HC';
    cPar.n=250;
    cPar.RigidConstraints=indRigid;
    [Vb_blob]=patchSmooth(F_blob,V_blob,[],cPar);
    indSmooth=unique(F_blob(:));
    V_blob(indSmooth,:)=Vb_blob(indSmooth,:);
    %Fix color data with new bottom surface
    Cb_blob=ones(size(Cb_blob));
    Cb_blob(all(ismember(Fb_blob,indRigid1),2))=2;

    meshStruct.nodes=V_blob;
end

%%
% Visualize blob mesh

hFig=cFigure; 
subplot(1,2,1); hold on;
gpatch(Fb_blob,V_blob,Cb_blob,'k',1);
% patchNormPlot(Fb_blob,V_blob);
% plotV(V_blob(indRigid,:),'g.','MarkerSize',25);
axisGeom(gca,fontSize);
colormap(gjet); icolorbar;
camlight headlight; 

hs=subplot(1,2,2); hold on; 
title('Cut view of solid mesh','FontSize',fontSize);
optionStruct.hFig=[hFig hs];
gpatch(Fb_blob,V_blob,'kw','none',0.25);
meshView(meshStruct,optionStruct);
axisGeom(gca,fontSize);
drawnow; 

%% Creating rigid body ground plate
% 

pointSpacingPlate=pointSpacingBlob;

% Compose outer curve of the plate
nPlateCurve=ceil((2*pi*plateRadius)/pointSpacingPlate);
t=linspace(0,2*pi,nPlateCurve);
t=t(1:end-1); 
x=plateRadius.*sin(t);
y=plateRadius.*cos(t); 
Vp_outer_curve=[x(:) y(:)];

% Create mesh out outer region
regionCell={Vp_outer_curve};
[F_plate,V_plate]=regionTriMesh2D(regionCell,pointSpacingPlate,0,0);
V_plate(:,3)=0; %Add z-direction

optionStruct.maxAngleDeviation=45;
[F_plate,V_plate]=tri2quadGroupSplit(F_plate,V_plate,optionStruct);

% Copy mesh for inner region
[F_plate_inner,V_plate_inner]=patchCleanUnused(fliplr(Fb_blob(Cb_blob==2,:)),V_blob);
[F_plate,V_plate,C_plate]=joinElementSets({F_plate,F_plate_inner},{V_plate,V_plate_inner});
[F_plate,V_plate]=mergeVertices(F_plate,V_plate);

center_of_mass_plate=mean(V_plate,1);

%%
% Visualizing plate mesh

cFigure; hold on;
gtitle('The plate surface mesh',fontSize);
gpatch(Fb_blob,V_blob,'kw','none',0.5);
gpatch(F_plate,V_plate,C_plate,'k',1);
patchNormPlot(F_plate,V_plate);
axisGeom(gca,fontSize);
colormap(gjet); icolorbar;
camlight headlight; 
drawnow; 

%% Creating rigid body shear surface

pointSpacingProbe=pointSpacingBlob/2; 

%Sketching side profile
x=[-probeLength-hemiSphereRadius -hemiSphereRadius -hemiSphereRadius];
y=[0 0 0];
z=[hemiSphereRadius*(1-probeOverlapFactor) hemiSphereRadius*(1-probeOverlapFactor) hemiSphereRadius*1.5];
V_probe_curve_sketch=[x(:) y(:) z(:)];

%Fillet sketch
np=100; %Number of points used to construct each fillet edge
[V_probe_curve]=filletCurve(V_probe_curve_sketch,filletProbe,np,0);
numPointsProbeCurve=ceil(max(pathLength(V_probe_curve))/pointSpacingProbe);
[V_probe_curve] = evenlySampleCurve(V_probe_curve,numPointsProbeCurve,'pchip',0);

% Extruding curve
% controlParametersExtrude.pointSpacing=pointSpacingProbe;
controlParametersExtrude.depth=hemiSphereRadius*2.5; 
controlParametersExtrude.numSteps=ceil(controlParametersExtrude.depth/pointSpacingProbe);
controlParametersExtrude.numSteps=controlParametersExtrude.numSteps+iseven(controlParametersExtrude.numSteps); %Force uneven
controlParametersExtrude.patchType='quad'; 
controlParametersExtrude.dir=0;
controlParametersExtrude.n=[0 1 0];
controlParametersExtrude.closeLoopOpt=0; 

[F_probe,V_probe]=polyExtrude(V_probe_curve,controlParametersExtrude);
F_probe=fliplr(F_probe); %Invert face orientation so normals point to blob

center_of_mass_probe=mean(V_probe,1);

%%
% Visualizing probe mesh

cFigure; hold on;
title('The probe surface mesh','fontSize',fontSize);
gpatch(Fb_blob,V_blob,'kw','none',0.5);
gpatch(F_plate,V_plate,'kw','none',0.5);
hl(1)=plotV(V_probe_curve_sketch,'k.-.','lineWidth',3,'MarkerSize',25);
hl(2)=plotV(V_probe_curve,'r-','lineWidth',3,'MarkerSize',25);
hl(3)=gpatch(F_probe,V_probe,'gw','k',1);
legend(hl,{'Sketched probe curve','Rounded probe curve','Probe surface mesh'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 
drawnow; 

%% Join model node sets

V=[V_blob; V_plate; V_probe];
F_plate=F_plate+size(V_blob,1);
F_probe=F_probe+size(V_blob,1)+size(V_plate,1);
Fb_all=[Fb_blob;F_plate;F_probe];

%%
% Visualizing model

cFigure; hold on;
gtitle('Model components',fontSize);
hl(1)=gpatch(Fb_blob,V,'rw','k',0.8);
hl(2)=gpatch(F_plate,V,'bw','k',0.8);
hl(3)=gpatch(F_probe,V,'gw','k',0.8);
legend(hl,{'Blob','Plate','Probe'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 
drawnow; 

%%

[N3,Vn]=patchNormal(Fb_blob,V);
[N1,N2]=vectorOrthogonalPair(N3); %Get orthogonal vector pair

%%
% Visualizing model

cFigure; hold on;
gtitle('Local material axes',fontSize);
hl(1)=gpatch(Fb_blob,V,'rw','k',0.9);
hl(2)=quiverVec(Vn,N1,pointSpacingBlob,'y');
hl(3)=quiverVec(Vn,N2,pointSpacingBlob,'g');
hl(4)=quiverVec(Vn,N3,pointSpacingBlob,'b');
legend(hl,{'Blob','1st axis','2nd axis','3rd axis'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 
drawnow; 

sdafasfs

%% Get contact surfaces
%

%%
% Visualize contact surfaces

cFigure; 
subplot(1,2,1); hold on;
title('Probe blob contact pair','fontsize',fontSize);
hl(1)=gpatch(F_probe,V,'rw','k',1);
patchNormPlot(F_probe,V);
hl(2)=gpatch(Fb_blob,V,'gw','k',1);
patchNormPlot(Fb_blob,V);
legend(hl,{'Master','Slave'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 

subplot(1,2,2); hold on;
title('Plate blob contact pair','fontsize',fontSize);
hl(1)=gpatch(F_plate,V,'rw','k',1);
patchNormPlot(F_plate,V);
hl(2)=gpatch(Fb_blob,V,'gw','k',1);
patchNormPlot(Fb_blob,V);
legend(hl,{'Master','Slave'}); clear hl;
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
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=step_size;
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
febio_spec.Material.material{2}.density=1;
febio_spec.Material.material{2}.center_of_mass=center_of_mass_plate;

febio_spec.Material.material{3}.ATTR.type='rigid body';
febio_spec.Material.material{3}.ATTR.id=3;
febio_spec.Material.material{3}.density=1;
febio_spec.Material.material{3}.center_of_mass=center_of_mass_probe;

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='hex8'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name='Blob'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(E_blob,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=E_blob;

febio_spec.Geometry.Elements{2}.ATTR.type='quad4'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat=2; %material index for this set 
febio_spec.Geometry.Elements{2}.ATTR.name='Plate'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id=size(E_blob,1)+(1:1:size(F_plate,1))'; %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL=F_plate;

febio_spec.Geometry.Elements{2}.ATTR.type='quad4'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat=3; %material index for this set 
febio_spec.Geometry.Elements{2}.ATTR.name='Probe'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id=size(E_blob,1)+size(F_plate,1)+(1:1:size(F_probe,1))'; %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL=F_probe;

% -> NodeSets


% % -> Surfaces
febio_spec.Geometry.Surface{1}.ATTR.name='contact_master1';
febio_spec.Geometry.Surface{1}.quad4.ATTR.lid=(1:1:size(F_plate,1))';
febio_spec.Geometry.Surface{1}.quad4.VAL=F_plate;

febio_spec.Geometry.Surface{2}.ATTR.name='contact_slave1';
febio_spec.Geometry.Surface{2}.quad4.ATTR.lid=(1:1:size(Fb_blob,1))';
febio_spec.Geometry.Surface{2}.quad4.VAL=Fb_blob;

febio_spec.Geometry.Surface{3}.ATTR.name='contact_master2';
febio_spec.Geometry.Surface{3}.quad4.ATTR.lid=(1:1:size(F_probe,1))';
febio_spec.Geometry.Surface{3}.quad4.VAL=F_probe;

febio_spec.Geometry.Surface{4}.ATTR.name='contact_slave2';
febio_spec.Geometry.Surface{4}.quad4.ATTR.lid=(1:1:size(Fb_blob(Cb_blob==1,:),1))';
febio_spec.Geometry.Surface{4}.quad4.VAL=Fb_blob(Cb_blob==1,:);

% -> Surface pairs
febio_spec.Geometry.SurfacePair{1}.ATTR.name='Contact1_plate_blob';
febio_spec.Geometry.SurfacePair{1}.master.ATTR.surface=febio_spec.Geometry.Surface{1}.ATTR.name;
febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;

febio_spec.Geometry.SurfacePair{2}.ATTR.name='Contact2_probe_blob';
febio_spec.Geometry.SurfacePair{2}.master.ATTR.surface=febio_spec.Geometry.Surface{3}.ATTR.name;
febio_spec.Geometry.SurfacePair{2}.slave.ATTR.surface=febio_spec.Geometry.Surface{4}.ATTR.name;

%Boundary condition section 
% -> Fix boundary conditions

% -> Prescribed boundary conditions on the rigid body
febio_spec.Boundary.rigid_body{1}.ATTR.mat=2;
febio_spec.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='x';
febio_spec.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='z';
febio_spec.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Rx';
febio_spec.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Ry';
febio_spec.Boundary.rigid_body{1}.fixed{6}.ATTR.bc='Rz';

febio_spec.Boundary.rigid_body{2}.ATTR.mat=3;
febio_spec.Boundary.rigid_body{2}.fixed{1}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{2}.fixed{2}.ATTR.bc='z';
febio_spec.Boundary.rigid_body{2}.fixed{3}.ATTR.bc='Rx';
febio_spec.Boundary.rigid_body{2}.fixed{4}.ATTR.bc='Ry';
febio_spec.Boundary.rigid_body{2}.fixed{5}.ATTR.bc='Rz';
febio_spec.Boundary.rigid_body{2}.prescribed.ATTR.bc='x';
febio_spec.Boundary.rigid_body{2}.prescribed.ATTR.lc=1;
febio_spec.Boundary.rigid_body{2}.prescribed.VAL=probeDisplacement;

%Contact section
% -> Contact 1
febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
febio_spec.Contact.contact{1}.ATTR.type='sticky';
febio_spec.Contact.contact{1}.penalty=10;
febio_spec.Contact.contact{1}.laugon=0;
febio_spec.Contact.contact{1}.tolerance=0.1;
febio_spec.Contact.contact{1}.minaug=1;
febio_spec.Contact.contact{1}.maxaug=10;
febio_spec.Contact.contact{1}.snap_tol=0;
febio_spec.Contact.contact{1}.max_traction=max_traction;
febio_spec.Contact.contact{1}.search_tolerance=0.1;

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
febioAnalysis.runMode='internal';%'internal';
febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=3; %Max log file checking time

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results 

if 1%runFlag==1 %i.e. a succesful run
    
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
    DN_magnitude=sqrt(sum(DN(:,3).^2,2));
    V_def=V+DN;
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
    X_DEF=V_DEF(:,1,:);
    Y_DEF=V_DEF(:,2,:);
    Z_DEF=V_DEF(:,3,:);
%     [CF]=vertexToFaceMeasure(Fb_all,DN_magnitude);
    
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; hold on;
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1=gpatch(Fb_blob,V_def,DN_magnitude,'k',1); %Add graphics object to animate
    hp2=gpatch(F_probe,V_def,'kw','none',0.5); %Add graphics object to animate
    hp3=gpatch(F_plate,V_def,'kw','none',0.5); %Add graphics object to animate
%     gpatch(Fb_all,V,0.5*ones(1,3),'none',0.25); %A static graphics object
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]); caxis manual;
    axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);
    camlight headlight;
    view(15,8); 
    drawnow; 
    % Set up animation features
    animStruct.Time=time_mat; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        DN=N_disp_mat(:,:,qt); %Current displacement
        DN_magnitude=sqrt(sum(DN.^2,2)); %Current displacement magnitude
        V_def=V+DN; %Current nodal coordinates
%         [CF]=vertexToFaceMeasure(Fb_all,DN_magnitude); %Current color data to use
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1 hp1 hp2 hp3]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,DN_magnitude,V_def,V_def}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    caxis([0 max(DN_magnitude)]); caxis manual;
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