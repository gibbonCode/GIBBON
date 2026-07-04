%% WIP_DEMO_gmsh_meshing
% Below is a demonstration for:
% 
% * Building geometry for a slab and dome using gmsh
% * Importing and visualizing the displacement results

%% Keywords
%
% * gmsh

%%
clear; close all; clc;

%%
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');
%path with pre-written .geo (gmsh) files
geo_path = fullfile(defaultFolder,'data','gmsh');
%list of prewritten .geo files used to generate meshes
geo_files = {'\rigid_dome.geo','\square_membrane.geo'};
%path to gmsh generated mesh files
mesh_path = savePath;
%list of gmsh generated mesh files
mesh_files = {'\rigid_dome.msh','\square_membrane.msh'};
%path to an installation of gmsh
gmsh_path = fullfile(defaultFolder,'lib_ext','gmsh','win64');

%%

parameter_names = {'membrane_L', 'membrane_W', 'membrane_H', 'radius'};
parameter_values = {50, 50, 0.5, 12.5};

%update the geo files 
for geo_ind = 1:length(geo_files)

    geo_text = txtfile2cell([geo_path,geo_files{geo_ind}]);
    
    for param_ind = 1:length(parameter_names)
        parameter_line =  find(startsWith(geo_text,parameter_names{param_ind}));
        if ~isempty(parameter_line)
            geo_text{parameter_line} = [parameter_names{param_ind}, ' = ', num2str(parameter_values{param_ind}),';'];
            
        end
    end
    
    cell2txtfile([savePath,geo_files{geo_ind}],geo_text);

end

%% Load Meshes
for i = 1:size(geo_files,2)
    runGmsh([savePath,geo_files{i}],gmsh_path)
end


%% Import Meshes
for i = 1:size(mesh_files,2)

Meshes{i} = read_gmsh([mesh_path, mesh_files{i}]);

end

%% Combine Meshes
MeshStruct = Meshes{1};
if size(Meshes,2)>1
    for i = 2:size(Meshes,2)

    vol_ct = max(MeshStruct.elementMaterialID);
    surf_ct = max(MeshStruct.boundaryMarker);
    
    MeshStruct_temp = Meshes{i};
    
    nodes = MeshStruct_temp.nodes;
    node_IDs = MeshStruct_temp.node_IDs+max(MeshStruct.node_IDs);
    facesBoundary = MeshStruct_temp.facesBoundary+max(MeshStruct.node_IDs);
    boundaryMarker = MeshStruct_temp.boundaryMarker+surf_ct;
    elements = MeshStruct_temp.elements+max(MeshStruct.node_IDs);
    elementMaterialID = MeshStruct_temp.elementMaterialID+vol_ct;
    part_name = MeshStruct_temp.loadNameStruct.MeshName;
    volumes_names = MeshStruct_temp.loadNameStruct.VolumeNames;
    facesBoundary_Names = MeshStruct_temp.loadNameStruct.SurfaceNames;
    

    MeshStruct.nodes = [MeshStruct.nodes;nodes];
    MeshStruct.node_IDs = [MeshStruct.node_IDs;node_IDs];
    MeshStruct.facesBoundary = [MeshStruct.facesBoundary;facesBoundary];
    MeshStruct.boundaryMarker = [MeshStruct.boundaryMarker;boundaryMarker];
    MeshStruct.elements = [MeshStruct.elements;elements];
    MeshStruct.elementMaterialID = [MeshStruct.elementMaterialID;elementMaterialID];
    MeshStruct.loadNameStruct.MeshName = {MeshStruct.loadNameStruct.MeshName;part_name};
    MeshStruct.loadNameStruct.VolumeNames = [MeshStruct.loadNameStruct.VolumeNames;volumes_names];
    MeshStruct.loadNameStruct.SurfaceNames = [MeshStruct.loadNameStruct.SurfaceNames;facesBoundary_Names];
    

    end
end

vol_ct = max(MeshStruct.elementMaterialID);
surf_ct = max(MeshStruct.boundaryMarker);


%% 
% % Plotting model boundary surfaces and a cut view
% 
E1=MeshStruct.elements; %The elements 
V1=MeshStruct.nodes; %The nodes (vertices)
Fb1=MeshStruct.facesBoundary; %The boundary faces
Cb1=MeshStruct.boundaryMarker; %The "colors" or labels for the boundary faces
%elementMaterialIndices=ones(size(E1,1),1); %Element material indices

%% Display the mesh

%Plot settings
fontSize=15;
faceAlpha1=1;
faceAlpha2=0.3;
markerSize=40;
markerSize2=20;
lineWidth=3;

hFig=cFigure;
hs1=subplot(1,2,1);
hold on;
for i = 1:surf_ct
    
    hl(i) = gpatch(Fb1(Cb1==i,:),V1,Cb1(Cb1==i),'k',faceAlpha1);

end
colormap(gjet(6))
legend(hl,MeshStruct.loadNameStruct.SurfaceNames,'Interpreter','none');
axisGeom(gca,fontSize);
drawnow;

hs2=subplot(1,2,2); hold on; 
title('Cut view of solid mesh','FontSize',fontSize);
optionStruct.hFig=[hFig hs2];
meshView(MeshStruct,optionStruct);
axisGeom(gca,fontSize);
C = findall(hFig,'type','ColorBar');
set(C,'TickLabels',MeshStruct.loadNameStruct.VolumeNames)
set(C,'TickLabelInterpreter','none')
%        'TickLabels',MeshStruct.loadNameStruct.VolumeNames,'TickLabelInterpreter','none')

drawnow;
%%

%%
% Defining file names
febioFebFileNamePart=['Gmsh_Example'];
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress
febioLogFileName_position=[febioFebFileNamePart,'_position_out.txt']; %Log file name for exporting stress
febioLogFileName_stress_prin=[febioFebFileNamePart,'_stress_prin_out.txt']; %Log file name for exporting principal stress
febioLogFileName_strain=[febioFebFileNamePart,'_strain_out.txt']; %Log file name for exporting stress
febioLogFileName_strain_prin=[febioFebFileNamePart,'_strain_prin_out.txt']; %Log file name for exporting principal stress


% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=30; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=25; %Maximum number of retires
dtmin=(1/numTimeSteps)/500; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size
runMode='external';

%Get a template with default settings 
[febio_spec]=febioStructTemplate;
febio_spec.ATTR.version='3.0'; 
febio_spec.Module.ATTR.type='solid'; 

%Control section
febio_spec.Control.analysis='STATIC';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.solver.max_ups=max_ups;
febio_spec.Control.solver.qnmethod='BFGS';
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;

% Create a default Material section
% Material parameter set
E_youngs1=0.1; %Material Young's modulus
nu1=0.48; %Material Poisson's ratio

defaultMaterial='defaultMaterial';
febio_spec.Material.material{1}.ATTR.name=defaultMaterial;
febio_spec.Material.material{1}.ATTR.type='rigid body';
febio_spec.Material.material{1}.ATTR.id=1;


%Output section 
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='sz';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{2}.ATTR.file=febioLogFileName_stress_prin;
febio_spec.Output.logfile.element_data{2}.ATTR.data='s1;s2;s3';
febio_spec.Output.logfile.element_data{2}.ATTR.delim=',';


%%

 febio_spec.Mesh.Nodes{1}.ATTR.name='all_nodes'; %The node set name
 febio_spec.Mesh.Nodes{1}.node.ATTR.id=MeshStruct.node_IDs; %The node id's
 febio_spec.Mesh.Nodes{1}.node.VAL=MeshStruct.nodes; %The nodel coordinates
 
 element_start = 1;
 for i = 1:vol_ct
    
    num_els = size(MeshStruct.elements(MeshStruct.elementMaterialID == i),1);
        
    febio_spec.Mesh.Elements{i}.ATTR.type = 'tet4'; %Element type
    febio_spec.Mesh.Elements{i}.ATTR.name = MeshStruct.loadNameStruct.VolumeNames{i}; %Name of this part
    febio_spec.Mesh.Elements{i}.elem.ATTR.id = [element_start:element_start+num_els-1]'; %Element id's
    febio_spec.Mesh.Elements{i}.elem.VAL = MeshStruct.elements(MeshStruct.elementMaterialID == i,:); %The element matrix
    febio_spec.MeshDomains.SolidDomain{i}.ATTR.name = MeshStruct.loadNameStruct.VolumeNames{i};
    febio_spec.MeshDomains.SolidDomain{i}.ATTR.mat = defaultMaterial;
    
    element_start = element_start+num_els;

 end
 
 
for i = 1:surf_ct
    
    element_ids = [1:size(MeshStruct.facesBoundary(MeshStruct.boundaryMarker == i,:),1)]';
    febio_spec.Mesh.Surface{i}.ATTR.name=MeshStruct.loadNameStruct.SurfaceNames{i}; %Name of this surface
    febio_spec.Mesh.Surface{i}.tri3.ATTR.id=element_ids; %Element id's
    febio_spec.Mesh.Surface{i}.tri3.VAL = MeshStruct.facesBoundary(MeshStruct.boundaryMarker == i,:); %The element matrix 

end

%%
matid = 1;
febio_spec.Material.material{matid}.ATTR.name='rigid';
febio_spec.Material.material{matid}.ATTR.type='rigid body';
febio_spec.Material.material{matid}.ATTR.id=matid;
febio_spec.Material.material{matid}.density = 1;
febio_spec.Material.material{matid}.center_of_mass = [0,0,0];

matid = 2;
febio_spec.Material.material{matid}.ATTR.name='membrane';
febio_spec.Material.material{matid}.ATTR.type='neo-Hookean';
febio_spec.Material.material{matid}.ATTR.id=matid;
febio_spec.Material.material{matid}.density = 1;
febio_spec.Material.material{matid}.E=100;
febio_spec.Material.material{matid}.v=.48;

febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat = 'rigid';
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat = 'membrane';
 %%
%% Constraints
%Rigid section
% ->Rigid body fix boundary conditions
febio_spec.Rigid.rigid_constraint{1}.ATTR.name='RigidFix';
febio_spec.Rigid.rigid_constraint{1}.ATTR.type='fix';
febio_spec.Rigid.rigid_constraint{1}.rb=1;
febio_spec.Rigid.rigid_constraint{1}.dofs='Rx,Ry,Rz,Ru,Rv,Rw';

%%Boundary conditions
% -> Fixed boundary conditions surfaces
bid = 1;
febio_spec.Boundary.bc{bid}.ATTR.name='fix_xyz_membrane';
febio_spec.Boundary.bc{bid}.ATTR.type='fix';
febio_spec.Boundary.bc{bid}.ATTR.node_set='@surface:membrane_fixed';
febio_spec.Boundary.bc{bid}.dofs='x,y,z';

bid = 2;
febio_spec.Boundary.bc{bid}.ATTR.name='fix_xyz_dome';
febio_spec.Boundary.bc{bid}.ATTR.type='fix';
febio_spec.Boundary.bc{bid}.ATTR.node_set='@surface:dome_fixed';
febio_spec.Boundary.bc{bid}.dofs='x,y,z';


%% Contacts
% -> Surface pairs
cid = 1;
febio_spec.Mesh.SurfacePair{cid}.ATTR.name= 'membrane_dome';
febio_spec.Mesh.SurfacePair{cid}.primary='membrane_bottom';
febio_spec.Mesh.SurfacePair{cid}.secondary='Dome';


febio_spec.Contact.contact{cid}.ATTR.type='sliding-facet-on-facet';
febio_spec.Contact.contact{cid}.ATTR.name= febio_spec.Mesh.SurfacePair{cid}.ATTR.name;
febio_spec.Contact.contact{cid}.ATTR.surface_pair=febio_spec.Mesh.SurfacePair{cid}.ATTR.name;
febio_spec.Contact.contact{cid}.two_pass=0;
febio_spec.Contact.contact{cid}.laugon= 1;
febio_spec.Contact.contact{cid}.tolerance=0.2;
febio_spec.Contact.contact{cid}.gaptol=0;
febio_spec.Contact.contact{cid}.minaug= 0;
febio_spec.Contact.contact{cid}.maxaug= 10;
febio_spec.Contact.contact{cid}.search_tol=0.01;
febio_spec.Contact.contact{cid}.search_radius=0;
%%febio_spec.Contact.contact{cid}.symmetric_stiffness=0;
febio_spec.Contact.contact{cid}.auto_penalty=0;
febio_spec.Contact.contact{cid}.penalty=1000;
febio_spec.Contact.contact{cid}.fric_coeff=0;



%% Loads

%LoadData section
% -> load_controller
lcid = 1;
febio_spec.LoadData.load_controller{lcid}.ATTR.id=lcid;
febio_spec.LoadData.load_controller{lcid}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{lcid}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{lcid}.points.point.VAL=[0 0; 0.2 0.05; 0.75 0.3; 1, 1];

febio_spec.Loads.surface_load{1}.ATTR.name='Membrane_Pressure';
febio_spec.Loads.surface_load{1}.ATTR.type='pressure';
febio_spec.Loads.surface_load{1}.ATTR.surface='membrane_top';
febio_spec.Loads.surface_load{1}.pressure.ATTR.lc=1;
febio_spec.Loads.surface_load{1}.pressure.VAL=0.67;
febio_spec.Loads.surface_load{1}.symmetric_stiffness=1;


%%
febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

%%
% 
febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; %Display information on the command window
febioAnalysis.runMode='internal';%'internal';

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%%
E1=MeshStruct.elements; %The elements 
V =MeshStruct.nodes; %The nodes (vertices)
Fb1=MeshStruct.facesBoundary; %The boundary faces
Cb1=MeshStruct.boundaryMarker; %The "colors" or labels for the boundary faces
elementMaterialIndices=MeshStruct.elementMaterialID; %Element material indices
%%
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
    % Importing element stress from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),1,1);     
    
    %Access data
    E_stress_mat=dataStruct.data;
    
            
    %%
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations
    
    [CV]=faceToVertexMeasure(E1,V,E_stress_mat(:,:,end));
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('$\sigma_{3}$ [MPa]','Interpreter','Latex')
    hp=gpatch(Fb1,V_DEF(:,:,end),CV,'k',1); %Add graphics object to animate
    %hp.Marker='.';
    hp.MarkerSize=markerSize2;
    hp.FaceColor='interp';
    
    %hp2=gpatch(E2,V_DEF(:,:,end),'w','none',0.5); %Add graphics object to animate
    
    axisGeom(gca,fontSize);
    colormap(flipud(gjet(250))); colorbar;
    caxis([min(E_stress_mat(:)) max(E_stress_mat(:))]);
    axis(axisLim(V_DEF)); %Set axis limits statically
    camlight headlight;
    
    % Set up animation features
    animStruct.Time=timeVec; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments
        
        [CV]=faceToVertexMeasure(E1,V,E_stress_mat(:,:,qt));
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),CV,V_DEF(:,:,qt)}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct); %Initiate animation feature
    drawnow;
        
end
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
