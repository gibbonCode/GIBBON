%% DEMO_febio_0006_sphere_indentation
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

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

%Specifying dimensions and number of elements for slab
sampleHeight=2; %Height
sampleWidth=sampleHeight; %Width 
sampleThickness=sampleHeight; %Thickness 
    
numElementsWidth=[4 3 5 4]; %Number of elemens in dir 1
numElementsThickness=numElementsWidth; %Number of elemens in dir 2
numElementsHeight=numElementsWidth; %Number of elemens in dir 3

%Material parameter set
youngsModuli=[0.3 10 0.3 10];
poissonsRatios=[0.4 0.1 0.4 0.1]; 

% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=10; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

%Contact parameters
contactInitialOffset=0.1;

boxOffsets=sampleHeight+contactInitialOffset;

%% Creating model geometry and mesh
% A box is created with tri-linear hexahedral (hex8) elements using the
% |hexMeshBox| function. The function offers the boundary faces with
% seperate labels for the top, bottom, left, right, front, and back sides.
% As such these can be used to define boundary conditions on the exterior. 

E=[];
elementMaterialID=[];
V=[];
Fb=[];
Cb=[];
for q=1:1:4
    
    % Create a box with hexahedral elements
    beamDimensions=[sampleWidth sampleThickness sampleHeight]; %Dimensions
    beamElementNumbers=[numElementsWidth(q) numElementsThickness(q) numElementsHeight(q)]; %Number of elements
    outputStructType=2; %A structure compatible with mesh view
    [meshStruct]=hexMeshBox(beamDimensions,beamElementNumbers,outputStructType);
    
    %Access elements, nodes, and faces from the structure
    E1=meshStruct.elements; %The elements
    V1=meshStruct.nodes; %The nodes (vertices)
    Fb1=meshStruct.facesBoundary; %The boundary faces
    Cb1=meshStruct.boundaryMarker; %The "colors" or labels for the boundary faces
    elementMaterialIndices=ones(size(E1,1),1); %Element material indices

    V1(:,3)=V1(:,3)+(q-1)*boxOffsets;
    
    E=[E;E1+size(V,1)];
    Fb=[Fb;Fb1+size(V,1)];
    V=[V;V1];    
    colorOffset=max(Cb(:));
    if isempty(colorOffset)
        colorOffset=0;
    end
    Cb=[Cb;Cb1+colorOffset];
    
    elementMaterialID=[elementMaterialID;q*ones(size(E1,1),1)];
    
end
V(:,3)=V(:,3)-min(V(:,3)); %Shift so bottom is at 0

%% 
% Plotting model boundary surfaces and a cut view

hFig=cFigure; 
hold on; 
title('Model boundary surfaces and labels','FontSize',fontSize);
gpatch(Fb,V,Cb,'k',faceAlpha1); 
colormap(gjet(250)); icolorbar;
axisGeom(gca,fontSize);
drawnow;

%%

logicTops=false(size(Cb,1),4);
logicBottoms=false(size(Cb,1),4);
for q=1:1:4
    logicTops(:,q)=Cb==6+(6*(q-1));   
    logicBottoms(:,q)=Cb==5+(6*(q-1));   
end


%% 
% Plotting model boundary surfaces 

hFig=cFigure; 
hold on; 
title('Contact faces','FontSize',fontSize);
gpatch(Fb,V,'kw','none',0.2);
for q=1:1:size(logicTops,2)
    gpatch(Fb(logicTops(:,q),:),V,Cb(logicTops(:,q),:),'k',1);
    gpatch(Fb(logicBottoms(:,q),:),V,Cb(logicBottoms(:,q),:),'k',1);
end

colormap(gjet(250)); icolorbar;
axisGeom(gca,fontSize);
drawnow;

%% Define boundary conditions

F_support=Fb(logicBottoms(:,1),:);

F_rigidBody=Fb(logicTops(:,end),:);
indNodesRigidBody=unique(F_rigidBody(:));
center_of_mass=mean(V(indNodesRigidBody,:),1);

%Supported nodes
bcSupportList=unique(F_support(:));

%Visualize BC's
hf=cFigure;
title('Boundary conditions model','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,V,'kw','none',faceAlpha2); 
hl(1)=gpatch(F_rigidBody,V,'rw','k',1); 
hl(2)=plotV(V(bcSupportList,:),'k.','MarkerSize',markerSize);
hl(3)=plotV(center_of_mass,'r.','MarkerSize',50);

legend(hl,{'Rigid body','BC support','Rigid body center of mass'});

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
for q=1:1:numel(youngsModuli)
    febio_spec.Material.material{q}.ATTR.type='neo-Hookean';
    febio_spec.Material.material{q}.ATTR.id=q;
    febio_spec.Material.material{q}.E=youngsModuli(q);
    febio_spec.Material.material{q}.v=poissonsRatios(q);
end

febio_spec.Material.material{numel(youngsModuli)+1}.ATTR.type='rigid body';
febio_spec.Material.material{numel(youngsModuli)+1}.ATTR.id=numel(youngsModuli)+1;
febio_spec.Material.material{numel(youngsModuli)+1}.density=1;
febio_spec.Material.material{numel(youngsModuli)+1}.center_of_mass=center_of_mass;

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
n=1;
for q=1:1:numel(youngsModuli)
    logicMaterialNow=(elementMaterialID==q);
    febio_spec.Geometry.Elements{q}.ATTR.type='hex8'; %Element type of this set
    febio_spec.Geometry.Elements{q}.ATTR.mat=q; %material index for this set
    febio_spec.Geometry.Elements{q}.ATTR.name=['Box_',num2str(q)]; %Name of the element set
    febio_spec.Geometry.Elements{q}.elem.ATTR.id=(n:1:(n-1+nnz(logicMaterialNow)))'; %Element id's
    febio_spec.Geometry.Elements{q}.elem.VAL=E(logicMaterialNow,:);    
    n=n+nnz(logicMaterialNow);
end

% -> NodeSets
febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcSupportList(:);

% -> Surfaces
c=1;
for q=2:1:numel(youngsModuli)
    F_contact_now=Fb(logicBottoms(q),:);
    febio_spec.Geometry.Surface{c}.ATTR.name=['contact_',num2str(c)];
    febio_spec.Geometry.Surface{c}.tri3.ATTR.lid=(1:1:size(F_contact_now,1))';
    febio_spec.Geometry.Surface{c}.tri3.VAL=F_contact_now;
    c=c+1;
end

for q=1:1:numel(youngsModuli)-1    
    c=numel(febio_spec.Geometry.Surface)+1;
    F_contact_now=Fb(logicTops(q),:);    
    febio_spec.Geometry.Surface{c}.ATTR.name=['contact_',num2str(c)];
    febio_spec.Geometry.Surface{c}.tri3.ATTR.lid=(1:1:size(F_contact_now,1))';
    febio_spec.Geometry.Surface{c}.tri3.VAL=F_contact_now;
end

% -> Surface pairs
surfaceInd=1;
for q=1:1:numel(youngsModuli)-1
    febio_spec.Geometry.SurfacePair{q}.ATTR.name=['Contact_',num2str(q)];
    febio_spec.Geometry.SurfacePair{q}.master.ATTR.surface=febio_spec.Geometry.Surface{q}.ATTR.name;
    febio_spec.Geometry.SurfacePair{q}.slave.ATTR.surface=febio_spec.Geometry.Surface{q+numel(youngsModuli)-1}.ATTR.name;
end

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.fix{1}.ATTR.bc='x';
febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{2}.ATTR.bc='y';
febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{3}.ATTR.bc='z';
febio_spec.Boundary.fix{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;

% -> Prescribed boundary conditions on the rigid body
febio_spec.Boundary.rigid_body{1}.ATTR.mat=numel(febio_spec.Material.material);
febio_spec.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='x';
febio_spec.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='Rx';
febio_spec.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Ry';
febio_spec.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Rz';
febio_spec.Boundary.rigid_body{1}.prescribed.ATTR.bc='z';
febio_spec.Boundary.rigid_body{1}.prescribed.ATTR.lc=1;
febio_spec.Boundary.rigid_body{1}.prescribed.VAL=bcPrescribeMagnitudes(3);

%Contact section
for q=1:1:numel(youngsModuli)-1
    febio_spec.Contact.contact{q}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{q}.ATTR.name;
    febio_spec.Contact.contact{q}.ATTR.type='sliding-elastic';
    febio_spec.Contact.contact{q}.penalty=1;
    febio_spec.Contact.contact{q}.auto_penalty=1;
    febio_spec.Contact.contact{q}.two_pass=1;
    febio_spec.Contact.contact{q}.laugon=1;
    febio_spec.Contact.contact{q}.tolerance=0.1;
    febio_spec.Contact.contact{q}.gaptol=0;
    febio_spec.Contact.contact{q}.minaug=0;
    febio_spec.Contact.contact{q}.maxaug=10;
    febio_spec.Contact.contact{q}.search_tol=0.01;
    febio_spec.Contact.contact{q}.search_radius=mean(pointSpacings)/2;
end

fsdfafdsfa

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
    DN_magnitude=sqrt(sum(DN(:,3).^2,2));
    V_def=V+DN;
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
    X_DEF=V_DEF(:,1,:);
    Y_DEF=V_DEF(:,2,:);
    Z_DEF=V_DEF(:,3,:);
    [CF]=vertexToFaceMeasure(Fb1,DN_magnitude);
    
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1=gpatch(Fb1,V_def,CF,'k',1); %Add graphics object to animate
    hp2=gpatch(E2,V_def,'kw','none',faceAlpha2); %Add graphics object to animate
    gpatch(Fb1,V,0.5*ones(1,3),'none',0.25); %A static graphics object
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]);
    axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);
    camlight headlight;
        
    % Set up animation features
    animStruct.Time=time_mat; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        DN=N_disp_mat(:,:,qt); %Current displacement
        DN_magnitude=sqrt(sum(DN.^2,2)); %Current displacement magnitude
        V_def=V+DN; %Current nodal coordinates
        [CF]=vertexToFaceMeasure(Fb1,DN_magnitude); %Current color data to use
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1 hp1 hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,CF,V_def}; %Property values for to set in order to animate
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
