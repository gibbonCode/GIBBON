%% DEMO_febio_0058_face_tube_loading
% Below is a demonstration for:
% 
% * Building triangulated surface geometry for a face
% * Meshing the face using pentahedral elements
% * Building model of a tube
% * Defining the boundary conditions 
% * Coding the febio structure
% * Running the model
% * Importing and visualizing results

%% Keywords
%
% * febio_spec version 2.5
% * febio, FEBio
% * face
% * contact, sliding, friction
% * pentahedral elements, penta6
% * static, solid
% * hyperelastic, Ogden
% * displacement logfile

%%

clear; close all; clc;

%%
% Plot settings
fontSize=10;
faceAlpha1=1;
faceAlpha2=0.3;
markerSize1=15;
markerSize2=10;
lineWidth=2; 

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
febioLogFileName_strainEnergy=[febioFebFileNamePart,'_energy_out.txt']; %Log file name for exporting strain energy density

volumeFactor=2; %Volume factor used in tetgen, larger means larger internal elements
layerThickness=5;
numSplit=1;

%Material parameter set
c1_1=1e-3; %Shear-modulus-like parameter
m1_1=2; %Material parameter setting degree of non-linearity
k_1=c1_1*100; %Bulk modulus

c1_2=c1_1*10; %Shear-modulus-like parameter
m1_2=2; %Material parameter setting degree of non-linearity
k_2=c1_2*100; %Bulk modulus

% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=35; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=10; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=(1/numTimeSteps); %Maximum time step size
symmetric_stiffness=0;
min_residual=1e-20;

runMode='internal';

%Boundary condition parameters
displacementMagnitude_y=0; %Displacement applied 
displacementMagnitude_z=-100; %Displacement applied 

%Contact parameters
contactPenalty=0.1;
laugon=0;
minaug=1;
maxaug=10;
fric_coeff=0.5;

%% Load face mesh

% [Fs,Vs]=graphicsModels(10);

[stlStruct] = import_STL(fullfile(defaultFolder,'data','STL','nefertiti_fine.stl'));
Fs=stlStruct.solidFaces{1}; %Faces
Vs=stlStruct.solidVertices{1}; %Vertices
[Fs,Vs]=mergeVertices(Fs,Vs); % Merging nodes


cFigure;  hold on; 
gpatch(Fs,Vs,'gw','k');
% gpatch(Eb,Vs,'bw','b',1,2);
% plotV(Vs(logicForce,:),'r.','MarkerSize',25);
axisGeom;
camlight headlight; 
drawnow; 

%%

Vs2=Vs;
Eb=patchBoundary(Fs,Vs);
cParSmooth.n=15;
cParSmooth.Method='LAP';
cParSmooth.RigidConstraints=unique(Eb(:));
Vs2=patchSmooth(Fs,Vs2,[],cParSmooth);

[~,~,N]=patchNormal(Fs,Vs2);
Vs2=Vs2-N*layerThickness; 
cParSmooth.n=5;
cParSmooth.Method='HC';
Vs2=patchSmooth(Fs,Vs2,[],cParSmooth);

V=[Vs;Vs2];
E_face=[Fs+size(Vs,1) Fs];

% [E_face,V]=subPenta(E_face,V,numSplit,3);

[F]=element2patch(E_face,V,'penta6');
Fp_tri=F{1};
Fp2=Fp_tri(1:size(Fs,1),:);
Fp1=Fp_tri(end-size(Fs,1)+1:end,:);

%%

cFigure; 
gpatch(Fs,Vs,'bw','none',1);
gpatch(F,V,'w','k',0.5);
axisGeom;
camlight headlight; 
drawnow; 

%%
cFigure; 
% gpatch(F,V,'kw','none',0.5);
% gpatch(Fp1,V,'rw','k',1);
gpatch(Fp1,V,'bw','k',1);
axisGeom;
camlight headlight; 
drawnow; 

%%

optionStruct.cylRadius=4;
optionStruct.numRadial=8;
optionStruct.cylHeight=200;
optionStruct.numHeight=[];
optionStruct.meshType='tri';
optionStruct.closeOpt=1;
[Fc,Vc,Cc]=patchcylinder(optionStruct);
R=euler2DCM([0 0.5*pi 0]);
Vc=Vc*R;
Vc(:,2)=Vc(:,2)-28;
Vc(:,3)=Vc(:,3)+20;

%%

[V_regions]=getInnerPoint(Fc,Vc); % Define region points
[regionA]=tetVolMeanEst(Fc,Vc); %Volume for regular tets

inputStruct.stringOpt='-pq1.2AaY';
inputStruct.Faces=Fc;
inputStruct.Nodes=Vc;
inputStruct.holePoints=[];
inputStruct.faceBoundaryMarker=Cc; %Face boundary markers
inputStruct.regionPoints=V_regions; %region points
inputStruct.regionA=regionA*volumeFactor;

% Mesh model using tetrahedral elements using tetGen 
[meshOutput]=runTetGen(inputStruct); %Run tetGen 

%% 
% Access model element and patch data
Fb_c=meshOutput.facesBoundary;
Cb_c=meshOutput.boundaryMarker;
Vc=meshOutput.nodes;
Ec=meshOutput.elements;

% Visualizing mesh using |meshView|, see also |anim8|
meshView(meshOutput);

%% Join node sets

Fb_c=Fb_c+size(V,1);
Ec=Ec+size(V,1);
V=[V;Vc];

%%
cFigure; 
gpatch(F,V,'rw','k',1);
gpatch(Fb_c,V,Cb_c,'k',1);
axisGeom;
colormap viridis;
icolorbar;
camlight headlight; 
drawnow; 

%%

%% Define contact surfaces

% The secondary surface
F_contact_secondary=fliplr(Fb_c(Cb_c==1,:));

Vc_secondary=patchCentre(F_contact_secondary,V);
Vc_Fp2=patchCentre(Fp1,V);
D=minDist(Vc_Fp2(:,[1 2]),Vc_secondary(:,[1 2]));

% The primary surface
F_contact_primary=Fp1(D<=3,:);

%%
% Visualize contact surfaces

cFigure; hold on;
title('Contact sets and normal directions','FontSize',fontSize);

gpatch(F,V,'kw','none',faceAlpha2); 

hl(1)=gpatch(F_contact_primary,V,'rw','k',1);
patchNormPlot(F_contact_primary,V);

hl(2)=gpatch(F_contact_secondary,V,'gw','k',1); 
patchNormPlot(F_contact_secondary,V);

legend(hl,{'Primary','Secondary'});

axisGeom(gca,fontSize);
camlight headlight;
drawnow;

%% Define boundary conditions

%Supported nodes
bcSupportList=unique(Fp2);

%Prescribed displacement nodes
bcPrescribeList=unique(Fb_c(Cb_c>1,:));

%%
% Visualize BC's

hf=cFigure; hold on;
title('Boundary conditions model','FontSize',fontSize);
gpatch(F,V,'kw','none',faceAlpha2); 
gpatch(Fb_c,V,'kw','none',faceAlpha2); 
hl2(1)=plotV(V(bcPrescribeList,:),'r.','MarkerSize',markerSize2);
hl2(2)=plotV(V(bcSupportList,:),'k.','MarkerSize',markerSize2);
legend(hl2,{'BC prescribe','BC support'});
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
febio_spec.Material.material{1}.c1=c1_1;
febio_spec.Material.material{1}.m1=m1_1;
febio_spec.Material.material{1}.c2=c1_1;
febio_spec.Material.material{1}.m2=-m1_1;
febio_spec.Material.material{1}.k=k_1;

febio_spec.Material.material{2}.ATTR.type='Ogden';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.c1=c1_2;
febio_spec.Material.material{2}.m1=m1_2;
febio_spec.Material.material{2}.c2=c1_2;
febio_spec.Material.material{2}.m2=-m1_2;
febio_spec.Material.material{2}.k=k_2;

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='penta6'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name='Face'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(E_face,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=E_face;

febio_spec.Geometry.Elements{2}.ATTR.type='tet4'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat=2; %material index for this set 
febio_spec.Geometry.Elements{2}.ATTR.name='Tube'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id=size(E_face,1)+(1:1:size(Ec,1))'; %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL=Ec;

% -> NodeSets
febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcSupportList(:);

febio_spec.Geometry.NodeSet{2}.ATTR.name='bcPrescribeList';
febio_spec.Geometry.NodeSet{2}.node.ATTR.id=bcPrescribeList(:);

% -> Surfaces
febio_spec.Geometry.Surface{1}.ATTR.name='contact_master';
febio_spec.Geometry.Surface{1}.tri3.ATTR.lid=(1:1:size(F_contact_secondary,1))';
febio_spec.Geometry.Surface{1}.tri3.VAL=F_contact_secondary;

febio_spec.Geometry.Surface{2}.ATTR.name='contact_slave';
febio_spec.Geometry.Surface{2}.tri3.ATTR.lid=(1:1:size(F_contact_primary,1))';
febio_spec.Geometry.Surface{2}.tri3.VAL=F_contact_primary;

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

% -> Prescribed boundary conditions
febio_spec.Boundary.prescribe{1}.ATTR.bc='y';
febio_spec.Boundary.prescribe{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Boundary.prescribe{1}.scale.ATTR.lc=1;
febio_spec.Boundary.prescribe{1}.scale.VAL=1;
febio_spec.Boundary.prescribe{1}.relative=1;
febio_spec.Boundary.prescribe{1}.value=displacementMagnitude_y;

febio_spec.Boundary.prescribe{2}.ATTR.bc='z';
febio_spec.Boundary.prescribe{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Boundary.prescribe{2}.scale.ATTR.lc=1;
febio_spec.Boundary.prescribe{2}.scale.VAL=1;
febio_spec.Boundary.prescribe{2}.relative=1;
febio_spec.Boundary.prescribe{2}.value=displacementMagnitude_z;

%Contact section
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
febio_spec.Contact.contact{1}.search_radius=0.1*sqrt(sum((max(V,[],1)-min(V,[],1)).^2,2)); 
febio_spec.Contact.contact{1}.symmetric_stiffness=0;
febio_spec.Contact.contact{1}.auto_penalty=1;
febio_spec.Contact.contact{1}.penalty=contactPenalty;
febio_spec.Contact.contact{1}.fric_coeff=fric_coeff;

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

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_strainEnergy;
febio_spec.Output.logfile.element_data{1}.ATTR.data='sed';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.element_data{1}.VAL=1:size(E_face,1);

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
febioAnalysis.runMode=runMode;%'internal';
febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

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
    
    V_def=V+DN;
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
    X_DEF=V_DEF(:,1,:);
    Y_DEF=V_DEF(:,2,:);
    Z_DEF=V_DEF(:,3,:);
  
    C=sqrt(sum(DN(:,3).^2,2));
        
    %%
    % Importing element strain energies from a log file
    [~,E_energy,~]=importFEBio_logfile(fullfile(savePath,febioLogFileName_strainEnergy)); %Element strain energy
    
    %Remove nodal index column
    E_energy=E_energy(:,2:end,:);
    
    %Add initial state i.e. zero energy
    sizImport=size(E_energy); 
    sizImport(3)=sizImport(3)+1;
    E_energy_mat_n=zeros(sizImport);
    E_energy_mat_n(:,:,2:end)=E_energy;
    E_energy=E_energy_mat_n;
    
    [FE_face,C_energy_face]=element2patch(E_face,E_energy(1:size(E_face,1),:,end),'penta6');
    [CV]=faceToVertexMeasure(FE_face,V,C_energy_face);
    
    %%
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1=gpatch(FE_face{1},V_def,CV,'none',1); %Add graphics object to animate
    hp1.FaceColor='Interp';
    hp2=gpatch(Fb_c,V_def,'k','none',0.3); %Add graphics object to animate
        
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([0 max(C_energy_face{1})/10]);
    axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);    
    camlight headlight;    
    lighting gouraud;
        
    % Set up animation features
    animStruct.Time=time_mat; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        DN=N_disp_mat(:,:,qt); %Current displacement        
        V_def=V+DN; %Current nodal coordinates
        
%         C=sqrt(sum(DN(:,3).^2,2)); %New color
        [FE_face,C_energy_face]=element2patch(E_face,E_energy(1:size(E_face,1),:,qt),'penta6');
        [CV]=faceToVertexMeasure(FE_face,V,C_energy_face);
                
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1 hp1 hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,CV,V_def}; %Property values for to set in order to animate
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
% Copyright (C) 2006-2022 Kevin Mattheus Moerman and the GIBBON contributors
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
