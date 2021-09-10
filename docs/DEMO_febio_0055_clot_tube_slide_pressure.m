%% DEMO_febio_0055_clot_tube_slide_pressure.m
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
% * febio_spec version 3.0
% * febio, FEBio
% * blood clot
% * contact, sliding, friction
% * rigid body constraints
% * hexahedral elements, hex8
% * quadrilaterl elements, quad4
% * shell elements
% * sphere
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
cMap=[1 0.5 0.4; 0.9 0.3 0.27; 0.8 0.2 0.18; 0.7 0.1 0.09; 0.6 0 0; 0.5 0 0; 0.4 0 0;];
[cMap]=resampleColormap(cMap,250);

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_strainEnergy=[febioFebFileNamePart,'_energy_out.txt']; %Log file name for exporting strain energy density

% Sphere parameters
sphereRadius=3;%
numElementsMantel=6;

% Ground plate parameters
tubeRadius=sphereRadius.*[1 0.1]; 
tubeAngle=3*(pi/180);
tubeLength=abs(diff(tubeRadius))/tan(tubeAngle);

% Material parameter set
c1=1e-4; %Shear-modulus-like parameter MPa
m1=2; %Material parameter setting degree of non-linearity
k_factor=10; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus
d=1e-9; %Density

% FEA control settings
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=15; %Optimum number of iterations
max_retries=25; %Maximum number of retires
symmetric_stiffness=0;
min_residual=1e-20;

timeTotal=1; %Analysis time
numTimeSteps=100; %Number of time steps desired
step_size=timeTotal/numTimeSteps;
dtmin=(timeTotal/numTimeSteps)/100; %Minimum time step size
dtmax=(timeTotal/numTimeSteps)*10; %Maximum time step size
analysisType='DYNAMIC';

%Contact parameters
contactPenalty=10;
laugon=0;
minaug=1;
maxaug=10;
fric_coeff=0.1; 

%Specifying load
pressureValue=7.5e-5;

%% Creating model geometry and mesh
% 

%Control settings
cPar.sphereRadius=sphereRadius;
cPar.coreRadius=cPar.sphereRadius/2;
cPar.numElementsMantel=numElementsMantel; 
cPar.numElementsCore=round(numElementsMantel*1.5); 
cPar.outputStructType=2;
cPar.makeHollow=0;
cPar.cParSmooth.n=25;

%Creating sphere
[meshOutput]=hexMeshSphere(cPar);

% Access model element and patch data
Fb_blob=meshOutput.facesBoundary;
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

pointSpacingBlob=mean(patchEdgeLengths(Fb_blob,V_blob));
pointSpacingTube=pointSpacingBlob;

rEnd=sphereRadius+(sphereRadius.*((sphereRadius-tubeRadius(2))/tubeLength));    
V_curve_tube=[sphereRadius rEnd 0; -tubeLength tubeRadius(2) 0;];

nResample=ceil(max(pathLength(V_curve_tube))./pointSpacingTube);
V_curve_tube=evenlySampleCurve(V_curve_tube,nResample,'pchip',0);

cPar.closeLoopOpt=1;
cPar.numSteps=[]; %If empty the number of steps is derived from point spacing of input curve
cPar.w=[1 0 0];
[F_tube,V_tube]=polyRevolve(V_curve_tube,cPar);
[F_tube,V_tube]=mergeVertices(F_tube,V_tube);
center_of_mass_tube=mean(V_tube,1);

%% Join model node sets

V=[V_blob; V_tube; ];
F_tube=F_tube+size(V_blob,1);

%%
% Visualizing model

cFigure; hold on;
gtitle('Model components',fontSize);
hl(1)=gpatch(Fb_blob,V,'rw','k',0.8);
hl(2)=gpatch(F_tube,V,'kw','k',0.5);
legend(hl,{'Blob','Tube'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 
drawnow; 

%% Get contact surfaces
%

F_contact_blob=Fb_blob;

%%
% Visualize contact surfaces

cFigure; hold on;
title('Tube blob contact pair','fontsize',fontSize);
hl(1)=gpatch(F_tube,V,'rw','k',0.8);
patchNormPlot(F_tube,V);
hl(2)=gpatch(F_contact_blob,V,'kw','k',0.5);
patchNormPlot(F_contact_blob,V);
legend(hl,{'Master','Slave'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 
drawnow; 

%% Get pressure surface

N=patchNormal(F_contact_blob,V);
x=[1 0 0];
D=dot(N,x(ones(size(N,1),1),:),2);

logicFace=D>-1e-6; %Logic for current face set
F_pressure=F_contact_blob(logicFace,:); %The current face set

%% 
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists. 

cFigure; hold on;
title('Pressure surface','fontsize',fontSize);

gpatch(F_pressure,V,'kw','k',0.5);
patchNormPlot(F_pressure,V);

% legend(hl,{'Master','Slave'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 
drawnow; 

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='3.0'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

%Control section
febio_spec.Control.analysis=analysisType;
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=timeTotal/numTimeSteps;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.solver.max_ups=max_ups;
febio_spec.Control.solver.symmetric_stiffness=symmetric_stiffness;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;

%Material section
materialName1='Material1';
febio_spec.Material.material{1}.ATTR.name=materialName1;
febio_spec.Material.material{1}.ATTR.type='Ogden unconstrained';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c1=c1;
febio_spec.Material.material{1}.m1=m1;
febio_spec.Material.material{1}.c2=c1;
febio_spec.Material.material{1}.m2=-m1;
febio_spec.Material.material{1}.cp=k;
febio_spec.Material.material{1}.density=d;

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
surfaceName1='Surface1_contact';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(F_tube,1))';
febio_spec.Mesh.Surface{1}.quad4.VAL=F_tube;

surfaceName2='Surface2_contact';
febio_spec.Mesh.Surface{2}.ATTR.name=surfaceName2;
febio_spec.Mesh.Surface{2}.quad4.ATTR.id=(1:1:size(F_contact_blob,1))';
febio_spec.Mesh.Surface{2}.quad4.VAL=F_contact_blob;

surfaceName3='Surface3_pressure';
febio_spec.Mesh.Surface{3}.ATTR.name=surfaceName3;
febio_spec.Mesh.Surface{3}.quad4.ATTR.id=(1:1:size(F_pressure,1))';
febio_spec.Mesh.Surface{3}.quad4.VAL=F_pressure;

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

%Loads
febio_spec.Loads.surface_load{1}.ATTR.type='pressure';
febio_spec.Loads.surface_load{1}.ATTR.surface=surfaceName3;
febio_spec.Loads.surface_load{1}.pressure.ATTR.lc=1;
febio_spec.Loads.surface_load{1}.pressure.VAL=pressureValue;
febio_spec.Loads.surface_load{1}.symmetric_stiffness=1;


%Rigid section 
% -> Prescribed rigid body boundary conditions
febio_spec.Rigid.rigid_constraint{1}.ATTR.name='RigidFix_1';
febio_spec.Rigid.rigid_constraint{1}.ATTR.type='fix';
febio_spec.Rigid.rigid_constraint{1}.rb=2;
febio_spec.Rigid.rigid_constraint{1}.dofs='Rx,Ry,Rz,Ru,Rv,Rw';

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
febio_spec.Contact.contact{1}.search_radius=0.01*sqrt(sum((max(V,[],1)-min(V,[],1)).^2,2)); 
febio_spec.Contact.contact{1}.symmetric_stiffness=0;
febio_spec.Contact.contact{1}.auto_penalty=1;
febio_spec.Contact.contact{1}.penalty=contactPenalty;
febio_spec.Contact.contact{1}.fric_coeff=fric_coeff;

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0; 1 1];

%Output section 
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_strainEnergy;
febio_spec.Output.logfile.element_data{1}.ATTR.data='sed';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.element_data{1}.VAL=1:size(E_blob,1);


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
febioAnalysis.runMode='external';%'internal';

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
    DN_magnitude=sqrt(sum(DN.^2,2));
    V_def=V+DN;
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
    X_DEF=V_DEF(:,1,:);
    Y_DEF=V_DEF(:,2,:);
    Z_DEF=V_DEF(:,3,:);
    %     [CF]=vertexToFaceMeasure(Fb_all,DN_magnitude);
    
    %%
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
    
    %%
    
    indBlob=unique(Fb_blob(:));
    t=linspace(0,2*pi,250)';
    
    V_def_blob=V(indBlob,:)+N_disp_mat(indBlob,:,end);
    
    [~,indMin]=min(V_def_blob(:,1));
    [~,indMax]=max(V_def_blob(:,1));
        
    xEnd=V_def_blob(indMin,1);        
    xStart=V_def_blob(indMax,1);        
    rEnd=sphereRadius+(xEnd.*((sphereRadius-tubeRadius(2))/tubeLength));    
    rStart=sphereRadius+(xStart.*((sphereRadius-tubeRadius(2))/tubeLength));    
    
    xMid=mean([xStart xEnd]);%sum([rStart rEnd].*[xStart xEnd])./sum([rStart rEnd]);
    rMid=sphereRadius+(xMid.*((sphereRadius-tubeRadius(2))/tubeLength));    
        
    V_plot_xEnd=[xEnd*ones(size(t)) rEnd*cos(t) rEnd*sin(t)];
    V_plot_xMid=[xMid*ones(size(t)) rMid*cos(t) rMid*sin(t)];
    V_plot_xStart=[xStart*ones(size(t)) rStart*cos(t) rStart*sin(t)];
       
   %%
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; hold on;
    ht=gtitle(['Radial stretch: ',num2str(rMid/sphereRadius)]);
    hp1=gpatch(Fb_blob,V_def,DN_magnitude,'none',1); %Add graphics object to animate
    
    hp2=plotV(V_plot_xEnd  ,'r-','LineWidth',3);    
    hp3=plotV(V_plot_xMid  ,'r-','LineWidth',3);    
    hp4=plotV(V_plot_xStart,'r-','LineWidth',3);    
    
    gpatch(F_tube,V_def,'kw','none',0.25); %Add graphics object to animate
    axisGeom(gca,fontSize); 
    colormap(cMap); colorbar;
    caxis([0 max(DN_magnitude(:))]); caxis manual;
    axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);
    camlight headlight; lighting gouraud;
    view(0,0);
%     view(-30,30); zoom(1.5);
    axis off;
    drawnow; 
    
    LMid=1;
    
    % Set up animation features
    animStruct.Time=time_mat; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        DN=N_disp_mat(:,:,qt); %Current displacement
        DN_magnitude=sqrt(sum(DN.^2,2));
        V_def=V+DN; %Current nodal coordinates
                
        V_def_blob=V_def(indBlob,:);
        [~,indMin]=min(V_def_blob(:,1));
        [~,indMax]=max(V_def_blob(:,1));
        
        xEnd=V_def_blob(indMin,1);
        xStart=V_def_blob(indMax,1);
        rEnd=sphereRadius+(xEnd.*((sphereRadius-tubeRadius(2))/tubeLength));
        rStart=sphereRadius+(xStart.*((sphereRadius-tubeRadius(2))/tubeLength));
        
        xMid=mean([xStart xEnd]);%sum([rStart rEnd].*[xStart xEnd])./sum([rStart rEnd]);
        rMid=sphereRadius+(xMid.*((sphereRadius-tubeRadius(2))/tubeLength));
        
        LMid=min(LMid,rMid/sphereRadius);
        
        V_plot_xEnd=[xEnd*ones(size(t)) rEnd*cos(t) rEnd*sin(t)];
        V_plot_xMid=[xMid*ones(size(t)) rMid*cos(t) rMid*sin(t)];
        V_plot_xStart=[xStart*ones(size(t)) rStart*cos(t) rStart*sin(t)];
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1 hp1 hp2 hp2 hp2 hp3 hp3 hp3 hp4 hp4 hp4 ht]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','XData','YData','ZData','XData','YData','ZData','XData','YData','ZData','String'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,DN_magnitude,...
            V_plot_xEnd(:,1),V_plot_xEnd(:,2),V_plot_xEnd(:,3),...
            V_plot_xMid(:,1),V_plot_xMid(:,2),V_plot_xMid(:,3),...
            V_plot_xStart(:,1),V_plot_xStart(:,2),V_plot_xStart(:,3),...
            ['Radial stretch: ',num2str(rMid/sphereRadius)]}; %Property values for to set in order to animate
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
% Copyright (C) 2006-2021 Kevin Mattheus Moerman and the GIBBON contributors
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
