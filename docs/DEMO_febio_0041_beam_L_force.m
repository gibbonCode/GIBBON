%% DEMO_febio_0041_beam_L_force
% Below is a demonstration for:
% 
% * Building geometry for a beam with hexahedral elements
% * Defining the boundary conditions 
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement results

%% Keywords
%
% * febio_spec version 4.0
% * febio, FEBio
% * beam force loading
% * force control boundary condition
% * hexahedral elements, hex8, hex20
% * beam, rectangular
% * static, solid
% * hyperelastic, Ogden
% * displacement logfile
% * stress logfile

%%

clear; close all; clc;

%%
% Plot settings
fontSize=15;
faceAlpha1=0.8;
faceAlpha2=1;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize=50;
markerSize2=20;
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
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

%Specifying dimensions and number of elements
L=85.4;
h=5.6;
b=50;
t=45;
phi=10;

sampleWidth=5.6;
sampleThickness=b; 
sampleHeight=L+h;

numElementsWidth=round(sampleWidth/h);
numElementsThickness=round(sampleThickness/h);
numElementsHeight=round(sampleHeight/h);

appliedForce=2800;

E_youngs1=6550; % 6.55e10 Pa = 65.5 GPa = 6550 MPa
v1=0.3;
E_youngs2=7200; %72 Gpa = 7200
v2=0.3;

nRefine=0;

% FEA control settings
numTimeSteps=20; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=10; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

runMode='external';% 'internal' or 'external'

%% CREATING MESHED BOX

%Create box 1
boxDim=[sampleWidth sampleThickness sampleHeight]; %Dimensions
boxEl=[numElementsWidth numElementsThickness numElementsHeight]; %Number of elements
[box1]=hexMeshBox(boxDim,boxEl);
E=box1.E;
V=box1.V;
Fb=box1.Fb;
Cb=box1.faceBoundaryMarker;

X=V(:,1); Y=V(:,2); Z=V(:,3);
VE=[mean(X(E),2) mean(Y(E),2) mean(Z(E),2)];

%%
% % Plotting boundary surfaces
% cFigure; hold on;
% title('Model surfaces','FontSize',fontSize);
% gpatch(Fb,V,Cb,'k',faceAlpha1);
% colormap(gjet(6)); icolorbar; 
% axisGeom;
% camlight headlight;
% set(gca,'FontSize',fontSize);
% drawnow; 

%% Make last element set "h" heigh

F_bottom=Fb(Cb==5,:);

logicBottomElements=any(ismember(E,F_bottom),2);
E_bottom=E(logicBottomElements,:);
[FE_bottom,~]=element2patch(E_bottom);

indV_FE_bottom=unique(FE_bottom(:));
mean_E_bottom=mean(V(indV_FE_bottom,:),1);

Z=V(:,3);
zMax=max(V(indV_FE_bottom,3));
zThreshold=zMax-(h/2);

indV_bottom=unique(F_bottom(:));
indV_bottomTop=indV_FE_bottom(V(indV_FE_bottom,3)>zThreshold);
% V(indV_bottomTop,3)=min(V(:,3))+h;

%% Find side faces to extrude

F_side=Fb(Cb==2,:);

logicSideBottom=all(ismember(F_side,indV_FE_bottom),2);
F_side_bottom=F_side(logicSideBottom,:);

layerThickness=t-h;
numStepsExtrude=ceil(layerThickness/h)+1;
dirSet=1; 
[Eq,Vq,Fq_start,Fq_end]=quadThick(F_side_bottom,V,dirSet,layerThickness,numStepsExtrude);

[Fq,~]=element2patch(Eq);

%% Merging element sets

F_top=Fb(Cb==6,:);
F_force=Fq_end;

F_force=F_force+size(V,1);
Eq=Eq+size(V,1);

elementMaterialIndices=[ones(size(E,1),1); 2*ones(size(Eq,1),1);];
E=[E;Eq];
V=[V;Vq];

[~,ind1,ind2]=unique(pround(V,5),'rows');
V=V(ind1,:);
E=ind2(E);
F_top=ind2(F_top);
F_force=ind2(F_force);

[FE,CE]=element2patch(E,elementMaterialIndices);

%%
% Plotting boundary surfaces
cFigure; hold on;
title('Model surfaces','FontSize',fontSize);

gpatch(FE,V,CE,'k',1);

gpatch(F_force,V,'y','y',1);
gpatch(F_top,V,'g','g',1);
% plotV(V,'k.','MarkerSize',25);

colormap(gjet(4)); icolorbar; 

axisGeom;
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow; 

%% Refining elements

C=[1:1:size(E,1)]';
if nRefine>0
    splitMethod=1;
    [E,V,C,CV]=subHex(E,V,nRefine,splitMethod);
end
elementMaterialIndices=elementMaterialIndices(C);

E1=E(elementMaterialIndices==1,:);
E2=E(elementMaterialIndices==2,:);

[FE,CF]=element2patch(E,C);
[D]=patchEdgeLengths(FE,V);

indTop=find(V(:,3)>(max(V(:,3))-max(D(:))/2));
indForce=find(V(:,1)>(max(V(:,1))-max(D(:))/2));


logicTopFaces=all(ismember(FE,indTop),2);
logicForceFaces=all(ismember(FE,indForce),2);

F_top=FE(logicTopFaces,:);
F_force=FE(logicForceFaces,:);
indBoundaryFaces=tesBoundary(FE,V);
Fb=FE(indBoundaryFaces,:);

%%
% Plotting boundary surfaces
cFigure; hold on;
title('Model surfaces','FontSize',fontSize);

gpatch(Fb,V,0.5*ones(1,3),'k',1);

gpatch(F_force,V,'r','k',1);
gpatch(F_top,V,'g','k',1);

plotV(V(indTop,:),'b.','MarkerSize',25);
plotV(V(indForce,:),'y.','MarkerSize',25);
% plotV(V,'k.','MarkerSize',25);

% colormap(gjet(6)); icolorbar; 
axisGeom;
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow; 

%% DEFINE BC's

%Supported nodes
bcSupportList=unique(F_top(:));

%Prescribed force nodes
bcPrescribeList=unique(F_force(:));
numForceNodes=numel(bcPrescribeList);

forceNormVec=[0 0 -1];
[R]=euler2DCM([0 phi/180*pi 0]);
forceNormVec=(R*forceNormVec')';

bcPrescribedForce=(appliedForce.*forceNormVec)/numForceNodes;

%%
% Visualize BC's

cFigure; hold on;
title('Boundary conditions','FontSize',fontSize);

gpatch(Fb,V,0.75*ones(1,3),'k',1);

plotV(V(bcSupportList,:),'b.','MarkerSize',markerSize);
plotV(V(bcPrescribeList,:),'r.','MarkerSize',markerSize);

quiverVec(V(bcPrescribeList,:),bcPrescribedForce(ones(numel(bcPrescribeList),1),:),25,'r');

axisGeom;
camlight headlight;
set(gca,'FontSize',fontSize);
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
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;

%Material section
materialName1='Material1';
febio_spec.Material.material{1}.ATTR.name=materialName1;
febio_spec.Material.material{1}.ATTR.type='neo-Hookean';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.E=E_youngs1;
febio_spec.Material.material{1}.v=v1;

materialName2='Material2';
febio_spec.Material.material{2}.ATTR.name=materialName2;
febio_spec.Material.material{2}.ATTR.type='neo-Hookean';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.E=E_youngs2;
febio_spec.Material.material{2}.v=v2;

%Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type 
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E1,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E1; %The element matrix

partName2='Part2';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='hex8'; %Element type 
febio_spec.Mesh.Elements{2}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=E2; %The element matrix

% -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcSupportList);

nodeSetName2='bcPrescribeList';
febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.VAL=mrow(bcPrescribeList);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat=materialName1;

febio_spec.MeshDomains.SolidDomain{2}.ATTR.name=partName2;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat=materialName2;

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_xyz';
febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.x_dof=1;
febio_spec.Boundary.bc{1}.y_dof=1;
febio_spec.Boundary.bc{1}.z_dof=1;

%Loads section
% -> Prescribed nodal forces
febio_spec.Loads.nodal_load{1}.ATTR.name='PrescribedForce';
febio_spec.Loads.nodal_load{1}.ATTR.type='nodal_force';
febio_spec.Loads.nodal_load{1}.ATTR.node_set=nodeSetName2;
febio_spec.Loads.nodal_load{1}.value.ATTR.lc=1;
febio_spec.Loads.nodal_load{1}.value.VAL=bcPrescribedForce;

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

febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';

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
