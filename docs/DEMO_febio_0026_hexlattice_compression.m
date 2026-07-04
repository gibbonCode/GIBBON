%% DEMO_febio_0026_hexlattice_compression
% Below is a demonstration for:
%
% * Building geometry for a cube with hexahedral elements
% * Defining the boundary conditions
% * Coding the febio structure
% * Running the model
% * Importing and visualizing the displacement and stress results

%% Keywords
%
% * febio_spec version 4.0
% * febio, FEBio
% * uniaxial loading
% * compression, tension, compressive, tensile
% * displacement control, displacement boundary condition
% * hexahedral elements, hex8
% * cube, box, rectangular
% * static, solid
% * hyperelastic, Ogden
% * displacement logfile
% * stress logfile

%%

clear; close all; clc;

%% Plot settings
% Plot settings
fontSize=15;
faceAlpha1=0.8;
faceAlpha2=1;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize=25;
cMap=gjet(4);

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
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress

%Specifying dimensions and number of elements
sampleSize=10;
shrinkFactor=0.2; 
nSplitPattern=1;

%Define applied displacement
appliedStrain=0.3; %Linear strain (Only used to compute applied stretch)
loadingOption='compression'; % or 'tension'
switch loadingOption
    case 'compression'
        stretchLoad=1-appliedStrain; %The applied stretch for uniaxial loading
    case 'tension'
        stretchLoad=1+appliedStrain; %The applied stretch for uniaxial loading
end
displacementMagnitude=(stretchLoad*sampleSize)-sampleSize; %The displacement magnitude

%Material parameter set
E_youngs1=0.1; %Material Young's modulus
nu1=0.4; %Material Poisson's ratio

% FEA control settings
numTimeSteps=10; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

runMode='external';% 'internal' or 'external'

%%
% Creating example geometry.
X=[-1;  1; 1; -1; -1;  1; 1; -1;];
Y=[-1; -1; 1;  1; -1; -1; 1;  1;];
Z=[-1; -1;-1; -1;  1;  1; 1;  1;];
V=sampleSize*[X(:) Y(:) Z(:)]/2;
E=1:8; %Element description of the 8-node cube (hexahedral element)

[E,V,C]=subHex(E,V,nSplitPattern); %Subdevide into 8 sub-cubes
[E,V,C]=hex2tet(E,V,C,1); %Convert to tetrahedral elements

[F,~]=element2patch(E,C); %Patch data for plotting
[indBoundary]=tesBoundary(F);

%%
% Create lattice structure

testCase=1;

controlParameter.meshType='hex'; %desired output mesh type
controlParameter.indBoundary=indBoundary; %indices of the boundary faces
switch testCase
    case 1
        controlParameter.shrinkFactor=shrinkFactor; %Strut sides are formed by shrinking the input mesh faces by this factor
        controlParameter.hexSplit=1; %Number of beam hex split steps
        controlParameter.latticeSide=2; %1=side 1 the edge lattice, 2=side 2 the dual lattice to the edge lattice
        [Es,Vs,Cs]=element2lattice(E,V,controlParameter);

    case 2

end
% Create patch Data for visualization
[Fs,CsF,CFs]=element2patch(Es,Cs); %Patch data for plotting

% Fs0=Fs(CsF==0,:);
% CsF0=CsF(CsF==0);

indB=tesBoundary(Fs);
Fb=Fs(indB,:);
CsFb=CsF(indB,:);
CFsb=CFs(indB,:);

%%
% Visualizing input mesh and lattic structures

cFigure;
hs=subplot(1,2,1);
title('The input mesh','fontSize',fontSize)
hold on;
gpatch(F,V,0.5*ones(1,3),'k',0.5);
axisGeom(gca,fontSize);
camlight headlight; lighting flat;

subplot(1,2,2);
title('Lattice side 1','fontSize',fontSize)
hold on;
gpatch(F,V,'none','k',0,2);
gpatch(Fb,Vs,'bw','k',0.5);
axisGeom(gca,fontSize);
camlight headlight; lighting flat;

drawnow;

%%

faceBoundaryMarker=zeros(size(Fb,1),1);

%Find top and bottom face sets
[Nb]=patchNormal(Fb,Vs);
zVec=[0 0 1];
d=dot(Nb,zVec(ones(size(Nb,1),1),:),2);
Z=Vs(:,3);
ZF=mean(Z(Fb),2);
logicTop_Fb=(d>0.9) & ZF>=(max(Vs(:,3))-eps(sampleSize));
logicBottom_Fb=(d<-0.9) & ZF<=(min(Vs(:,3))+eps(sampleSize));

xVec=[1 0 0];
d=dot(Nb,xVec(ones(size(Nb,1),1),:),2);
X=Vs(:,1);
XF=mean(X(Fb),2);
logicSides_Fb1=(d>0.9) & XF>=(max(Vs(:,1))-eps(sampleSize));
logicSides_Fb2=(d<-0.9) & XF<=(min(Vs(:,1))+eps(sampleSize));

yVec=[0 1 0];
d=dot(Nb,yVec(ones(size(Nb,1),1),:),2);
Y=Vs(:,2);
YF=mean(Y(Fb),2);
logicSides_Fb3=(d>0.9) & YF>=(max(Vs(:,2))-eps(sampleSize));
logicSides_Fb4=(d<-0.9) & YF<=(min(Vs(:,2))+eps(sampleSize));

faceBoundaryMarker(logicBottom_Fb)=1;
faceBoundaryMarker(logicTop_Fb)=2;
faceBoundaryMarker(logicSides_Fb1)=3;
faceBoundaryMarker(logicSides_Fb2)=4;
faceBoundaryMarker(logicSides_Fb3)=5;
faceBoundaryMarker(logicSides_Fb4)=6;

%% Smoothen lattice

%Get "clean" surface mesh for only the boundary
[Fb_clean,Vb_clean,indFix]=patchCleanUnused(Fb,Vs);

%Find indices of points to hold on to during smoothing
indKeep=Fb_clean(faceBoundaryMarker>0,:);
indKeep=unique(indKeep(:));

%Smoothing
cPar.Method='HC';
cPar.n=15;
cPar.RigidConstraints=indKeep;
[Vb_clean]=tesSmooth(Fb_clean,Vb_clean,[],cPar);

%Override coordinates in full mesh with smoothed coordinates
ind=Fb(:);
ind=unique(ind(:));
Vs(ind,:)=Vb_clean;

%%
% Visualize smoothed mesh

cFigure; hold on;
title('Smoothed mesh','FontSize',fontSize);
gpatch(F,V,'kw','none',0.1);
gpatch(Fb,Vs,faceBoundaryMarker,'k',1);
% plotV(Vb_clean(indKeep,:),'r.','markerSize',15)
axisGeom;
camlight headlight;
set(gca,'FontSize',fontSize);
% axis off
drawnow;

%% DEFINE BC's

%Supported nodes
bcSupportList=unique(Fb(faceBoundaryMarker==1,:));

%Prescribed force nodes
bcPrescribeList=unique(Fb(faceBoundaryMarker==2,:));

%%
% Visualize BC's

cFigure; hold on;
title('Boundary conditions','FontSize',fontSize);
gpatch(Fb,Vs,'kw','none',0.4);
hl(1)=plotV(Vs(bcSupportList,:),'k.','MarkerSize',markerSize);
hl(2)=plotV(Vs(bcPrescribeList,:),'r.','MarkerSize',markerSize);
legend(hl,{'BC full support','BC prescribed Z displacement'})
axisGeom;
camlight headlight;
set(gca,'FontSize',fontSize);
drawnow;

%%

logicNotBeam=CsF==1;
Fss=Fs(logicNotBeam,:);

controlParGroup.outputType='label';
g=tesgroup(Fss,controlParGroup);

R=nan(max(g),1);

for q=1:1:max(g)
    logicNow= (g==q);
    F_now=Fss(logicNow,:);
    indNow=unique(F_now(:));
    indNow=indNow(ismember(indNow,Fb));    
    V_now=Vs(indNow,:);
    V_now=V_now-mean(V_now);
    R(q)=mean(sqrt(sum(V_now.^2,2)));
end
strutThicknessMean=2*mean(R);

%%
% cFigure; hold on;
% gpatch(F,V,'none','k',0,2);
% gpatch(Fb,Vs,'none','k',1,1);
% gpatch(Fss,Vs,g,'k',0.5);
% axisGeom(gca,fontSize);
% camlight headlight; lighting flat;
% % colormap spectral; icolorbar
% drawnow;

%% Check porosity

vol_lattice=sum(hexVol(Es,Vs)); %Volume of hexahedra
porosity_lattice=vol_lattice./sampleSize.^3; %Porosity

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
febio_spec.Material.material{1}.v=nu1;

%Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(Vs,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=Vs; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(Es,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=Es; %The element matrix

% -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcSupportList);

nodeSetName2='bcPrescribeList';
febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.VAL=mrow(bcPrescribeList);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_xyz';
febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.x_dof=1;
febio_spec.Boundary.bc{1}.y_dof=1;
febio_spec.Boundary.bc{1}.z_dof=1;

febio_spec.Boundary.bc{2}.ATTR.name='zero_displacement_xy';
febio_spec.Boundary.bc{2}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{2}.x_dof=1;
febio_spec.Boundary.bc{2}.y_dof=1;
febio_spec.Boundary.bc{2}.z_dof=0;

febio_spec.Boundary.bc{3}.ATTR.name='prescibed_displacement_z';
febio_spec.Boundary.bc{3}.ATTR.type='prescribed displacement';
febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{3}.dof='z';
febio_spec.Boundary.bc{3}.value.ATTR.lc=1;
febio_spec.Boundary.bc{3}.value.VAL=displacementMagnitude;
febio_spec.Boundary.bc{3}.relative=0;

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

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='s1';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

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
    Vs_DEF=N_disp_mat+repmat(Vs,[1 1 size(N_disp_mat,3)]);

    %%
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations

    DN_magnitude=sqrt(sum(N_disp_mat(:,:,end).^2,2)); %Current displacement magnitude

    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('Displacement magnitude [mm]','Interpreter','Latex')
    hp=gpatch(Fb,Vs_DEF(:,:,end),DN_magnitude,'k',1); %Add graphics object to animate

    hp.FaceColor='interp';

    axisGeom(gca,fontSize);
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]);
    axis(axisLim(Vs_DEF)); %Set axis limits statically
    camlight headlight;

    % Set up animation features
    animStruct.Time=timeVec; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments
        DN_magnitude=sqrt(sum(N_disp_mat(:,:,qt).^2,2)); %Current displacement magnitude

        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={Vs_DEF(:,:,qt),DN_magnitude}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct); %Initiate animation feature
    drawnow;

    %%
    % Importing element stress from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),0,1);

    %Access data
    E_stress_mat=dataStruct.data;
    E_stress_mat(isnan(E_stress_mat))=0;

    %%
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations

    [CV]=faceToVertexMeasure(Es,Vs,E_stress_mat(:,:,end));

    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('$\sigma_{zz}$ [MPa]','Interpreter','Latex')
    hp=gpatch(Fb,Vs_DEF(:,:,end),CV,'k',1); %Add graphics object to animate

    hp.FaceColor='interp';

    axisGeom(gca,fontSize);
    colormap(gjet(250)); colorbar;
    caxis([min(E_stress_mat(:)) max(E_stress_mat(:))]/3);
    axis(axisLim(Vs_DEF)); %Set axis limits statically
    camlight headlight;

    % Set up animation features
    animStruct.Time=timeVec; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments

        [CV]=faceToVertexMeasure(Es,Vs,E_stress_mat(:,:,qt));

        %Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
        animStruct.Set{qt}={Vs_DEF(:,:,qt),CV}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct); %Initiate animation feature
    drawnow;

    %%
    % Importing nodal forces from a log file

    [dataStruct]=importFEBio_logfile(fullfile(savePath,febioLogFileName_force),0,1); %Nodal forces

    %Access data
    timeVec=dataStruct.time;
    f_sum_x=squeeze(sum(dataStruct.data(bcPrescribeList,1,:),1));
    f_sum_y=squeeze(sum(dataStruct.data(bcPrescribeList,2,:),1));
    f_sum_z=squeeze(sum(dataStruct.data(bcPrescribeList,3,:),1));

    %%
    % Visualize force data

    displacementApplied=timeVec.*displacementMagnitude;

    cFigure; hold on;
    xlabel('$u$ [mm]','Interpreter','Latex');
    ylabel('$F_z$ [N]','Interpreter','Latex');
    hp=plot(displacementApplied(:),f_sum_z(:),'b-','LineWidth',3);
    grid on; box on; axis square; axis tight;
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
