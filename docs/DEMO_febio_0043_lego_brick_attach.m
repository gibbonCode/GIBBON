%% DEMO_febio_0043_lego_brick_attach.m
% Below is a demonstration for:
% 
% * Building geometry for a lego block with tetrahedral elements% 
% * Defining the boundary conditions 
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the results

%% Keywords
%
% * febio_spec version 2.5
% * febio, FEBio% * 
% * contact, sliding, sticky, friction
% * tetrahedral elements, tet4
% * triangular elements, tri3
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
surfaceModelPath=fullfile(defaultFolder,'data','libSurf');
surfaceModelName='lego_block_2x4_035.mat';

volumeFactor=500;

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement

% Material parameter set
c1=1e-3; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
k_factor=25; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus

% FEA control settings
numTimeSteps=25; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=12; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/500; %Minimum time step size
dtmax=(1/numTimeSteps); %Maximum time step size
symmetric_stiffness=0;
min_residual=1e-20;

% Contact parameters
contactPenalty=0.1;
laugon=0;
minaug=1;
maxaug=10;
fric_coeff=0.1;

% Boundary condition parameters
displacementMagnitude=-1;
initialOffset=0.1;

%% Creating model geometry and mesh

% Import surface model
surfaceData=load(fullfile(surfaceModelPath,surfaceModelName));

%Access model data
F_brick=surfaceData.F;
V_brick=surfaceData.V;
C_brick=surfaceData.C;

% Merging nodes
[F_brick,V_brick]=mergeVertices(F_brick,V_brick);

% Smooth mesh
logicBoundaryNodes=false(size(V_brick,1),1);
for q=unique(C_brick)'
    Eb=patchBoundary(F_brick(C_brick==q,:),V_brick);    
    logicBoundaryNodes(Eb(:))=1;
end

cPar.Method='HC';
cPar.n=10;
cPar.RigidConstraints=find(logicBoundaryNodes);
[V_brick]=patchSmooth(F_brick,V_brick,[],cPar); 

% Get mesh point spacing
pointSpacing=mean(patchEdgeLengths(F_brick,V_brick));

% Rotate model
R=euler2DCM([-0.5*pi 0 0]);
V_brick=V_brick*R;

% Center model
V_brick=V_brick-mean(V_brick,1); %Center on mean
V_brick(:,3)=V_brick(:,3)-min(V_brick(:,3)); %Shift so bottom is at 0 in z-dir

%% 
% Plotting model boundary surfaces and a cut view

cFigure; hold on;
title('Model boundary surfaces and labels','FontSize',fontSize);
gpatch(F_brick,V_brick,C_brick,'k',faceAlpha1); 
plotV(V_brick(logicBoundaryNodes,:),'r.','MarkerSize',25);
colormap(gjet(250)); icolorbar;
axisGeom(gca,fontSize);
camlight headlight; 
drawnow;

%%

colorIndicesTop=[1:3 32:35 37];
colorIndicesBottom=[4:11];

cFigure; 
title('Model boundary surfaces and labels','FontSize',fontSize);
gpatch(F_brick,V_brick,C_brick,'none',0.25); 

logicPlot1=ismember(C_brick,colorIndicesBottom);
gpatch(F_brick(logicPlot1,:),V_brick,C_brick(logicPlot1,:),'k',1); 

colormap(gjet(250)); icolorbar;
axisGeom(gca,fontSize);
camlight headlight; 

drawnow;

cFigure; 
title('Model boundary surfaces and labels','FontSize',fontSize);
gpatch(F_brick,V_brick,C_brick,'none',0.25); 

logicPlot1=ismember(C_brick,colorIndicesTop);
gpatch(F_brick(logicPlot1,:),V_brick,C_brick(logicPlot1,:),'k',1); 

colormap(gjet(250)); icolorbar;
axisGeom(gca,fontSize);
camlight headlight; 

drawnow;

%% Mesh using tetrahedral elements

stringOpt='-pq1.2AaY';

inputStruct.stringOpt=stringOpt;
inputStruct.Faces=F_brick;
inputStruct.Nodes=V_brick;
inputStruct.holePoints=[];
inputStruct.faceBoundaryMarker=C_brick; %Face boundary markers
inputStruct.regionPoints=getInnerPoint(F_brick,V_brick); %region points
inputStruct.regionA=tetVolMeanEst(F_brick,V_brick)*volumeFactor; %Volume for regular tets
inputStruct.minRegionMarker=2; %Minimum region marker
 
% Mesh model using tetrahedral elements using tetGen 
[meshOutput]=runTetGen(inputStruct); %Run tetGen 

% Access model element and patch data
Fb=meshOutput.facesBoundary;
Cb=meshOutput.boundaryMarker;
V=meshOutput.nodes;
CE=meshOutput.elementMaterialID;
E=meshOutput.elements;

%% Visualizing mesh using |meshView|, see also |anim8|dsfaasdfs
meshView(meshOutput);

%%

Fb1=Fb;
Fb2=Fb;
V1=V; 
V2=V;
V2(:,3)=V2(:,3)+max(V2(:,3))+initialOffset;

%%
% Plotting joined geometry
cFigure;
title('Joined node sets','FontSize',fontSize);
hold on;
gpatch(Fb1,V1,Cb,'k',faceAlpha1); 
gpatch(Fb2,V2,Cb,'k',faceAlpha1); 

colormap(gjet(6)); icolorbar; 
axisGeom(gca,fontSize);
camlight headlight;
drawnow;

%% Joining node and element sets
V=[V1;V2;]; %Combined node sets
Fb2=Fb2+size(V1,1);
E1=E;
E2=E+size(V1,1); 
E=[E1;E2]; 

%%
% Plotting joined geometry
cFigure;
title('Joined node sets','FontSize',fontSize);
hold on;
gpatch(Fb1,V,'r','k'); 
gpatch(Fb2,V,'b','k');
axisGeom(gca,fontSize);
camlight headlight;
drawnow;

%% Define contact surfaces

% The rigid master surface of the sphere
F_contact_master=fliplr(Fb2(ismember(Cb,colorIndicesTop),:));

% The deformable slave surface of the slab
F_contact_slave=fliplr(Fb1(ismember(Cb,colorIndicesBottom),:));

% Plotting surface models
cFigure; hold on;
title('Contact sets and normal directions','FontSize',fontSize);

gpatch(Fb1,V,'kw','none',0.5); 
gpatch(Fb2,V,'kw','none',0.5); 

hl(1)=gpatch(F_contact_master,V,'gw','k',1); 
patchNormPlot(F_contact_master,V);
hl(2)=gpatch(F_contact_slave,V,'bw','k',1);
patchNormPlot(F_contact_slave,V);

legend(hl,{'Master','Slave'});

axisGeom(gca,fontSize);
camlight headlight;
drawnow;

%% Define boundary conditions

% Supported nodes
F_base=Fb1(ismember(Cb,37),:);
bcSupportList=unique(F_base);

% Prescribed displacement nodes
F_top=Fb2(ismember(Cb,12),:);
bcPrescribeList=unique(F_top);

%%
% Visualize BC's
hf=cFigure;
title('Boundary conditions model','FontSize',fontSize);
hold on;

gpatch(Fb1,V,'kw','none',faceAlpha2); 
gpatch(Fb2,V,'kw','none',faceAlpha2); 
gpatch(F_base,V,'rw','k',1); 
gpatch(F_top,V,'gw','k',1); 

clear hl;
hl(1)=plotV(V(bcSupportList,:),'r.','MarkerSize',markerSize);
hl(2)=plotV(V(bcPrescribeList,:),'g.','MarkerSize',markerSize);
legend(hl,{'BC support','BC prescribe'});

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

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='tet4'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name='Propeller'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=E;

% -> NodeSets
febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcSupportList(:);

febio_spec.Geometry.NodeSet{2}.ATTR.name='bcPrescribeList';
febio_spec.Geometry.NodeSet{2}.node.ATTR.id=bcPrescribeList(:);

% -> Surfaces
febio_spec.Geometry.Surface{1}.ATTR.name='contact_master';
febio_spec.Geometry.Surface{1}.tri3.ATTR.lid=(1:1:size(F_contact_master,1))';
febio_spec.Geometry.Surface{1}.tri3.VAL=F_contact_master;

febio_spec.Geometry.Surface{2}.ATTR.name='contact_slave';
febio_spec.Geometry.Surface{2}.tri3.ATTR.lid=(1:1:size(F_contact_slave,1))';
febio_spec.Geometry.Surface{2}.tri3.VAL=F_contact_slave;

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

febio_spec.Boundary.fix{4}.ATTR.bc='x';
febio_spec.Boundary.fix{4}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Boundary.fix{5}.ATTR.bc='y';
febio_spec.Boundary.fix{5}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;

% -> Prescribe boundary conditions
febio_spec.Boundary.prescribe{1}.ATTR.bc='z';
febio_spec.Boundary.prescribe{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Boundary.prescribe{1}.scale.ATTR.lc=1;
febio_spec.Boundary.prescribe{1}.scale.VAL=1;
febio_spec.Boundary.prescribe{1}.relative=1;
febio_spec.Boundary.prescribe{1}.value=displacementMagnitude;

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
febioAnalysis.runMode='internal';%'internal';
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
        
    % Importing element strain energies from a log file
    [~,E_energy,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_strainEnergy)); %Element stresses
    
    %Remove nodal index column
    E_energy=E_energy(:,2:end,:);
    
    %Add initial state i.e. zero displacement
    sizImport=size(E_energy); 
    sizImport(3)=sizImport(3)+1;
    E_energy_mat_n=zeros(sizImport);
    E_energy_mat_n(:,:,2:end)=E_energy;
    E_energy=E_energy_mat_n;
        
    C=faceToVertexMeasure(E,V,E_energy(:,:,end));
        
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    % Create basic view and store graphics handle to initiate animation    
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1=gpatch(Fb,V_def,C,'k',1); %Add graphics object to animate    
    hp1.FaceColor='interp';
    hp2=gpatch(Fc,V_def,'kw','none',0.5); %Add graphics object to animate
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([0 max(E_energy(:))/40]);
    axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);
    camlight headlight;
    
    % Set up animation features
    animStruct.Time=time_mat; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        DN=N_disp_mat(:,:,qt); %Current displacement
        DN_magnitude=sqrt(sum(DN.^2,2)); %Current displacement magnitude
        V_def=V+DN; %Current nodal coordinates
        
        C=faceToVertexMeasure(E,V,E_energy(:,:,qt));
         
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1 hp1 hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,C,V_def}; %Property values for to set in order to animate
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
% Copyright (C) 2019  Kevin Mattheus Moerman
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
