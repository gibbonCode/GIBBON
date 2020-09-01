%% DEMO_febio_0027_layer_spatially_varying_material
% Below is a demonstration for:
% 
% * Building geometry for a cube with hexahedral elements
% * Defining the boundary conditions 
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement and stress results

%% Keywords
%
% * febio_spec version 2.5
% * febio, FEBio
% * uniaxial loading
% * compression, tension, compressive, tensile
% * displacement control, displacement boundary condition
% * hexahedral elements, hex8
% * cube, box, rectangular
% * static, solid
% * hyperelastic, Ogden
% * displacement logfile
% * SED logfile

%%

clear; close all; clc;

%% Plot settings
% Plot settings
fontSize=15;
faceAlpha1=0.8;
faceAlpha2=1;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize1=10;
cMap=gjet(4);

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_sed=[febioFebFileNamePart,'_sed_out.txt']; %Log file name for exporting strain energy density

% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

%% DEFINING AND VISUALIZING THE PARAMETER MAP

numElemIncLayer=3;
numElemTopLayer=3;
numElemLayers=numElemTopLayer+numElemIncLayer;
displacementMagnitude=-0.3.*numElemLayers;

[G]=textImage('GIBBON','Arial',25,5);
G=flipud(G);
G=G-min(G(:));
G=G./max(G(:));

S=zeros(size(G,1),size(G,2),numElemLayers);
S(:,:,1:numElemIncLayer)=repmat(G,[1 1 numElemIncLayer]);

%%
% Control parameters

nBins=50; 
minC=1e-3; %minimum value
maxC=minC*100; %Maximum value
c1_range_ini=linspace(minC,maxC,nBins); %Value range
k_factor=50;

%% 
% VISUALIZING THE MAPPING

[F,V,C]=ind2patch(true(size(S)),S,'vb'); 
[C_rgb]=gray2RGBColorMap(C,jet(250),[min(S(:)) max(S(:))]);

[Fs1,Vs1,Cs1]=ind2patch(S>0,S,'vb'); 
[Fs2,Vs2,Cs2]=ind2patch(S==0,S,'vb'); 

cFigure;

subplot(1,2,1);
title('Stiff inclusion','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
patch('Faces',Fs1,'Vertices',Vs1,'FaceColor','flat','CData',Cs1,'EdgeColor','k','FaceAlpha',1);
axis equal; view(3); axis tight; axis vis3d; grid on; view([-20,22]);
colormap(cMap); caxis([min(S(:)) max(S(:))]); colorbar; 
set(gca,'FontSize',fontSize);

subplot(1,2,2);
title('Soft matrix','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize); hold on;
patch('Faces',Fs2,'Vertices',Vs2,'FaceColor','flat','CData',Cs2,'EdgeColor','k','FaceAlpha',1);
axis equal; view(3); axis tight; axis vis3d; grid on; view([-20,22]);
colormap(cMap); caxis([min(S(:)) max(S(:))]); colorbar;
set(gca,'FontSize',fontSize);

drawnow;

%% BUILD MODEL

%Create hexahedral elements with function based colors
[E,V,C]=ind2patch(true(size(S)),S,'hu'); 

%Define element parameter mapping
elementMaterialID=C;
elementMaterialID=elementMaterialID-min(elementMaterialID(:));
elementMaterialID=elementMaterialID./max(elementMaterialID(:)); %Normalized
elementMaterialID=round(elementMaterialID.*(nBins-1))+1; %1-nPar

indUni=unique(elementMaterialID(:)); %Unique indices of used materials
c1=c1_range_ini(indUni); %Select relevant points
numMaterials=numel(c1);

%Fix indices 
indFix1=1:numel(indUni);
indFix2=zeros(nBins,1);
indFix2(indUni)=indFix1;
elementMaterialID=indFix2(elementMaterialID);

%Reorder elementMaterialIndices and element matrix
[elementMaterialID,indSort]=sort(elementMaterialID);
E=E(indSort,:);

[F,PF]=element2patch(E,elementMaterialID);

%Get boundary faces for light plotting
[indBoundary]=tesBoundary(F,V);
Fb=F(indBoundary,:);

%% SET UP BOUNDARY CONDITIONS

%List of nodes for applying displacement
logicTopNodes=abs(V(:,3)-max(V(:,3)))<=max(eps(V(:,3)));
bcPrescribeList=find(logicTopNodes);

%List of nodes to fix
logicBottomNodes=abs(V(:,3)-min(V(:,3)))<=max(eps(V(:,3)));
bcRigidList=find(logicBottomNodes);

%%
% Visualize BC's

cFigure; hold on;
title('Boundary conditions','FontSize',fontSize);
gpatch(Fb,V,'kw','none',0.4);
hl(1)=plotV(V(bcRigidList,:),'k.','MarkerSize',markerSize1);
hl(2)=plotV(V(bcPrescribeList,:),'r.','MarkerSize',markerSize1);
legend(hl,{'BC full support','BC prescribed pressure'})
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
febio_spec.ATTR.version='2.5'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

%Control section
febio_spec.Control.analysis.ATTR.type='static';
febio_spec.Control.title='Lattice analysis';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax; 
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;
febio_spec.Control.max_refs=max_refs;
febio_spec.Control.max_ups=max_ups;

%Material section
for q=1:1:numMaterials
    febio_spec.Material.material{q}.ATTR.type='Ogden';
    febio_spec.Material.material{q}.ATTR.id=q;
    febio_spec.Material.material{q}.c1=c1(q);
    febio_spec.Material.material{q}.m1=2;
    febio_spec.Material.material{q}.c2=c1(q);
    febio_spec.Material.material{q}.m2=-2;
    febio_spec.Material.material{q}.k=k_factor.*c1(q);
end

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
n=1;
for q=1:1:numMaterials
    logicMaterialNow=(elementMaterialID==q);
    febio_spec.Geometry.Elements{q}.ATTR.type='hex8'; %Element type of this set
    febio_spec.Geometry.Elements{q}.ATTR.mat=q; %material index for this set
    febio_spec.Geometry.Elements{q}.ATTR.name=['Layer_mat',num2str(q)]; %Name of the element set
    febio_spec.Geometry.Elements{q}.elem.ATTR.id=(n:1:(n-1+nnz(logicMaterialNow)))'; %Element id's
    febio_spec.Geometry.Elements{q}.elem.VAL=E(logicMaterialNow,:);    
    n=n+nnz(logicMaterialNow);
end

% -> NodeSets
febio_spec.Geometry.NodeSet{1}.ATTR.name='bcRigidList';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcRigidList(:);

febio_spec.Geometry.NodeSet{2}.ATTR.name='bcPrescribeList';
febio_spec.Geometry.NodeSet{2}.node.ATTR.id=bcPrescribeList(:);

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.fix{1}.ATTR.bc='x';
febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{2}.ATTR.bc='y';
febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{3}.ATTR.bc='z';
febio_spec.Boundary.fix{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
% febio_spec.Boundary.fix{4}.ATTR.bc='x';
% febio_spec.Boundary.fix{4}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
% febio_spec.Boundary.fix{5}.ATTR.bc='y';
% febio_spec.Boundary.fix{5}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;

% -> Prescribe boundary conditions
febio_spec.Boundary.prescribe{1}.ATTR.bc='z';
febio_spec.Boundary.prescribe{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{2}.ATTR.name;
febio_spec.Boundary.prescribe{1}.scale.ATTR.lc=1;
febio_spec.Boundary.prescribe{1}.scale.VAL=1;
febio_spec.Boundary.prescribe{1}.relative=1;
febio_spec.Boundary.prescribe{1}.value=displacementMagnitude;

%Output section 
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_sed;
febio_spec.Output.logfile.element_data{1}.ATTR.data='sed';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.element_data{1}.VAL=1:size(E,1);

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
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_sed),1,1);
    
    %Access data
    E_energy=dataStruct.data;
    
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    [CV]=faceToVertexMeasure(E,V,E_energy(:,:,end));
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('$\sigma_{zz}$ [MPa]','Interpreter','Latex')
    hp=gpatch(Fb,V_DEF(:,:,end),CV,'k',1); %Add graphics object to animate

    hp.FaceColor='interp';
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([min(E_energy(:)) max(E_energy(:))]/8);    
    axis(axisLim(V_DEF)); %Set axis limits statically    
    camlight headlight;        
        
    % Set up animation features
    animStruct.Time=timeVec; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        
        [CV]=faceToVertexMeasure(E,V,E_energy(:,:,qt));
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;
    
    %% 
    % Calculate metrics to visualize stretch-stress curve
    
    DZ_set=N_disp_mat(bcPrescribeList,end,:); %Z displacements of the prescribed set
    DZ_set=mean(DZ_set,1); %Calculate mean Z displacements across nodes
    stretch_sim=(DZ_set(:)+sampleHeight)./sampleHeight; %Derive stretch
    stress_cauchy_sim=mean(squeeze(E_energy(:,end,:)),1)';
    
    %%    
    % Visualize stress-stretch curve
    
    cFigure; hold on;    
    title('Uniaxial stress-stretch curve','FontSize',fontSize);
    xlabel('$\lambda$ [.]','FontSize',fontSize,'Interpreter','Latex'); 
    ylabel('$\sigma_{zz}$ [MPa]','FontSize',fontSize,'Interpreter','Latex'); 
    
    plot(stretch_sim(:),stress_cauchy_sim(:),'r-','lineWidth',lineWidth);
    
    view(2); axis tight;  grid on; axis square; box on; 
    set(gca,'FontSize',fontSize);
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
