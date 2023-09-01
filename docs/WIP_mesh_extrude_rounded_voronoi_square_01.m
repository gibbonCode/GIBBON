clear; close all; clc;

%%
% Plot Settings
fontSize=25;
edgeWidth=3;
markerSize=35;
lineWidthVec=2;

%% Control parameters

w=1; % Width of square space
cellSpacing=0.2; % Cell spacing
pointSpacing=cellSpacing/15; % Mesh point spacing
layerThickness=pointSpacing*10; %Layer thickness for the mesh
pointSpacingThickness=pointSpacing*3; %Point spacing in the thickness direction

smoothWeight=0.5; %0 =linear, 1=fully smoothed (spline)
groupingOptionStruct.outputType='label';
scaleFactor=0.85;

shearAngle=0;
a=1./cosd(shearAngle); %Scale x, e.g. to compensate for shear induced thickness reduction
b=1; %Scale y
S=[a 0 0; ...
   tand(-shearAngle) b      0;...
   0 0      1];

randomOpt=1;

uniNoiseWidth=cellSpacing/3; 
stdNoise=cellSpacing/3; 

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath'))); 
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress sigma_z

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

%% Derived parameters
numThickness=ceil((layerThickness/pointSpacingThickness));

Si=inv(S);

%% Create

V_square=w*[0 0 0; 1 0 0; 1 1 0; 0 1 0]; %Square; 
V_square=V_square*S;

[V]=evenlySpaceCurve(V_square,cellSpacing,'linear',1,1:1:size(V_square,1));

[~,indCorners_V]=minDist(V_square,V);

if randomOpt==1
    %Top
    indSelect=find(V(:,2)>max(V(:,2)-eps(1)));
    indSelect=indSelect(~ismember(indSelect,indCorners_V));
    V(indSelect,1)=V(indSelect,1)+uniNoiseWidth*(rand(numel(indSelect),1)-0.5);
    
    %Bottom
    indSelect=find(V(:,2)<min(V(:,2)+eps(1)));
    indSelect=indSelect(~ismember(indSelect,indCorners_V));
    V(indSelect,1)=V(indSelect,1)+uniNoiseWidth*(rand(numel(indSelect),1)-0.5);
    
    %Right
    indSelect=find(V(:,1)>max(V(:,1)-eps(1)));
    indSelect=indSelect(~ismember(indSelect,indCorners_V));
    V(indSelect,2)=V(indSelect,2)+uniNoiseWidth*(rand(numel(indSelect),1)-0.5);
    
    %Left
    indSelect=find(V(:,1)<min(V(:,1)+eps(1)));
    indSelect=indSelect(~ismember(indSelect,indCorners_V));
    V(indSelect,2)=V(indSelect,2)+uniNoiseWidth*(rand(numel(indSelect),1)-0.5);
    
    %%
    
%     cFigure; hold on;
%     title('Edge points'); 
%     plotV(V,'k.','markerSize',50,'LineWidth',3);
%     plotV(V(indCorners_V,:),'r.','markerSize',60,'LineWidth',3);
%     plotV(V(indSelect,:),'g.','markerSize',75,'LineWidth',3);
%     
%     axisGeom; view(2);
%     drawnow;
    
end
%%

%Create Delaunay derived mesh
regionCell={V(:,[1 2])};

if randomOpt==1
    [Fs,Vs]=regionTriMeshRand2D(regionCell,cellSpacing,stdNoise*ones(1,2),0,0);
else
    [Fs,Vs]=regionTriMesh2D(regionCell,cellSpacing,0,0);
end

Es=patchBoundary(Fs);
indBoundary_s=unique(Es(:));

[Vd,Fd,~,Cd]=patch_dual(Vs,Fs,1);            

Vd(:,3)=0;
Vs(:,3)=0;

%%

cFigure; hold on;
title('Deformed triangulation and cell mesh');
hp1=gpatch(Fs,Vs,'w','k',1,edgeWidth);
hp2=gpatch(Fd,Vd,'bw','b',0.25,edgeWidth);
% legend([hp1 hp2],{'Input triangulation','Cell mesh'},'Location','NorthOutside');
axis equal tight;
set(gca,'FontSize',fontSize)
drawnow;

%%

Vd_def=Vd;
Vs_def=Vs;
Vs=Vs*Si;
Vd=Vd*Si;
V_square=V_square*Si; 

%%

Fds=Fd;
Vds=[];
for q=1:1:numel(Fd)
    f=Fd{q};
    c=Cd{q};
    v=Vs(c,:);
    
    [fs,vs]=scalePatch(f,Vd,scaleFactor,v);
    
    Fds{q}=fs+size(Vds,1);
    Vds=[Vds;vs];
end

%%

cFigure; hold on;
title('Undeformed geometry'); 
hp1=gpatch(Fs,Vs,'w','k',1,edgeWidth);
hp2=gpatch(Fds,Vds,'bw','b',0.5,edgeWidth);
% legend([hp1 hp2],{'Triangulation','Shrunk cell mesh'},'Location','NorthOutside');
axis equal tight;
set(gca,'FontSize',fontSize)
drawnow;

%%

[Fq,Vq,Fc,Vc]=dualClad(Fs,Vs,1-scaleFactor,1);
[Fbt,Vbt,Cbt]=joinElementSets({[Fq(:,[1 2 3]);Fq(:,[3 4 1])],Fc},{Vq,Vc});
[Fbt,Vbt]=mergeVertices(Fbt,Vbt);

V_total=[Vds;Vbt];
Fbt=Fbt+size(Vds,1);

%%
Fdt=[];
for q=1:1:numel(Fds)
    f=Fds{q};
    v=patchCentre(f,V_total);
    e=patchEdges(f,0);
    ind=repmat(1:1:size(f,1),size(f,2),1);
    f=[e size(V_total,1)+ind(:)];
    Fdt=[Fdt;f];
    V_total=[V_total;v];
end

F_total=[Fbt;Fdt;];
C_total=[ones(size(Fbt,1),1); 2*ones(size(Fdt,1),1);];

[F_total,V_total]=mergeVertices(F_total,V_total);

%%

E_branches=patchBoundary(F_total(C_total==1,:));
[G_branches,~,groupSize]=tesgroup(E_branches,groupingOptionStruct);
[~,indMax]=max(groupSize);

E_branchBoundary=E_branches(G_branches==indMax,:);
[~,indCorners]=minDist(V_square,V_total);


E_cells=patchBoundary(F_total(C_total==2,:));
[G2]=tesgroup(E_cells,groupingOptionStruct);

%%

Eb=patchBoundary(F_total);
logicMember=ismember(E_branchBoundary,Eb(:));
logicSideEdges=all(logicMember,2);
E_branchBoundary_sides=E_branchBoundary(logicSideEdges,:);
E_branchBoundary_inner=E_branchBoundary(~logicSideEdges,:);

logicMember=ismember(E_cells,Eb(:));
logicSideEdges=all(logicMember,2);
E_cell_boundary=E_cells(logicSideEdges,:);

logicMember=ismember(E_cells,E_branchBoundary(:));
logicSideEdges=~any(logicMember,2);
E_cell_inner=E_cells(logicSideEdges,:);

%%

E1=E_branchBoundary_sides;
G1=(1:1:size(E1,1))';
maxIndexBranchSides=max(G1(:));

E2=E_branchBoundary_inner;
G2=tesgroup(E2,groupingOptionStruct)+max(G1(:));
maxIndexBranches=max(G2(:));

E3=E_cell_boundary;
G3=(1:1:size(E3,1))'+max(G2(:));
maxIndexCellBoundary=max(G3(:));

E4=E_cell_inner;
G4=tesgroup(E4,groupingOptionStruct)+max(G3(:));

E_total=[E1;E2;E3;E4];
G_total=[G1;G2;G3;G4];

%%

cFigure; 
subplot(1,2,1); hold on;
title('Coarse cell and ECM triangulation');
gpatch(F_total(C_total==1,:),V_total,'rw','k',1,1);
gpatch(F_total(C_total==2,:),V_total,'bw','k',1,1);
axis equal tight;
set(gca,'FontSize',fontSize)

subplot(1,2,2); hold on;
title('Boundary edge sets');
GV=faceToVertexMeasure(E_total,V_total,G_total);
gpatch(E_total,V_total,GV,'interp',1,3);
% plotV(VC(indCorners,:),'k.','MarkerSize',markerSize);
colormap(gjet); %icolorbar;
axis equal tight;
set(gca,'FontSize',fontSize)

drawnow;

%%

E_raw=E_total;
V_raw=V_total;

%% Process resampling
for q=1:1:max(G_total(:))
    
    logicNow=G_total==q;
    E_now=E_total(logicNow,:); %Get current edges to resample
    
    E_total=E_total(~logicNow,:); %Remove edges from set
    G_total=G_total(~logicNow); %Remove group from set
    
    [E_new,V_new]=resampleEdgeSet(E_now,V_total,pointSpacing,smoothWeight);
    
    V_total=[V_total; V_new]; %Add new points
    E_total=[E_total; E_new];
    G_total=[G_total; q*ones(size(E_new,1),1)];
    
end

%% Remove unused points and force merging of resampled contours

[E_total,V_total]=patchCleanUnused(E_total,V_total);


%% Remove unused points and force merging of resampled contours

[E_total,V_total]=patchCleanUnused(E_total,V_total);
[E_total,V_total]=mergeVertices(E_total,V_total);

%%

cFigure; hold on;
title('Smoothed boundary sets');
GV=faceToVertexMeasure(E_total,V_total,G_total);
gpatch(E_raw,V_raw,'none','k',1,1);
gpatch(E_total,V_total,GV,'interp',1,edgeWidth);
% plotV(VC,'k.','MarkerSize',markerSize)

axis tight; axis equal; colormap(gjet); %icolorbar;
set(gca,'FontSize',fontSize);
axis off
gdrawnow;

%% Create logic arrays to capture the edges for different boundary groupings relevant for meshing

logicBranches=ismember(G_total,1:maxIndexBranches);
logicBoundaryCells=ismember(G_total,maxIndexBranchSides+1:maxIndexCellBoundary);
logicInnerCells=G_total>maxIndexCellBoundary;

%%

cFigure; hold on;

hp1=gpatch(E_total(logicBranches,:),V_total,'none','b',1,4);
hp2=gpatch(E_total(logicBoundaryCells,:),V_total,'none','r',1,3);
hp3=gpatch(E_total(logicInnerCells,:),V_total,'none','g',1,2);
legend([hp1 hp2 hp3],{'Branch boundary','Boundary cells','Interior cells'});
axis tight; axis equal;
set(gca,'FontSize',fontSize);
gdrawnow;

%% Set up mesh regions
% The single closed outer boundary is defined as a region with the cells as
% "holes" in this boundary. This defines the connective tissue mesh region.
% Next the inner cells are defined as mesh regions followed by the outer
% boundary cells.

indNow=edgeListToCurve(E_total(logicBranches,:));
indNow=indNow(1:end-1);
Es=E_total(logicInnerCells,:);
Gs=tesgroup(Es,groupingOptionStruct);
R=cell(1,max(Gs)+1);
R{1,1}=V_total(indNow,1:2);
regionSpec={};
for q=1:1:max(Gs)
    E_now=Es(Gs==q,:);
    indNow=edgeListToCurve(E_now);
    indNow=indNow(1:end-1);
    R{1,q+1}=V_total(indNow,1:2);
    regionSpec{q+1}={V_total(indNow,1:2);};
end
regionSpec{1}=R;

Es=E_total(logicBoundaryCells,:);
Gs=tesgroup(Es,groupingOptionStruct);
for q=1:1:max(Gs)
    E_now=Es(Gs==q,:);
    indNow=edgeListToCurve(E_now);
    indNow=indNow(1:end-1);
    regionSpec{end+1}={V_total(indNow,1:2);};
end

%% Using multiRegionTriMesh2D to create the triangulated mesh

[F,V,regionLabelSurface]=multiRegionTriMesh2D(regionSpec,pointSpacing,0,0);
V(:,3)=0; %Add z coordinate

Eb=patchBoundary(F); %The outer boundary edges
Ebb=patchBoundary(F(regionLabelSurface==1,:)); %The connective tissue boundary edges

%%

cMap=[0.5 0.5 0.5; gjet(max(regionLabelSurface)-1)]; %A custom colormap that starts grey so connective tissue stands out

cFigure; hold on;
gpatch(F,V,regionLabelSurface,'k',1,0.5);
% gpatch(E_raw,V_raw,'none','k',1,4);
gpatch(Eb,V,'none','k',1,4);
gpatch(Ebb,V,'none','k',1,3);
axis tight; axis equal; colormap(cMap);
% icolorbar;
set(gca,'FontSize',fontSize);
axis off;
drawnow;

%% Define fibre directions for the mesh

F_con=F(regionLabelSurface==1,:);
V_F_con=patchCentre(F_con,V);

logicCon=~all(ismember(Ebb,Eb),2);
V_Ebb=patchCentre(Ebb(logicCon,:),V);
N_Ebb=vecnormalize(V(Ebb(logicCon,1),:)-V(Ebb(logicCon,2),:));

[~,indMin]=minDist(V_F_con,V_Ebb);
N_con1=N_Ebb(indMin,:);

% nz=[0 0 1];
% N_con2=cross(N_con1,nz(ones(size(N_con1,1),1),:),2);

%%

cFigure; hold on;
gpatch(F,V,regionLabelSurface,'none',0.5);
% gpatch(Eb,V,'none','k',1,4);
gpatch(Ebb,V,'none','k',1,3);
quiverLine(V_F_con,N_con1,pointSpacing,'b',lineWidthVec); % plotV(V_F_con,'k.');
% quiverLine(V_F_con,N_con2,pointSpacing,'r',lineWidthVec); % plotV(V_F_con,'k.');
axis tight; axis equal; colormap(cMap);
set(gca,'FontSize',fontSize);
drawnow;

%% Thicken the mesh to 3D volumetric elements
[Ep,Vp,Fp1,Fp2]=patchThick(F,V,1,layerThickness,numThickness);

%Copy region label so it holds for solid mesh
regionLabelElements=repmat(regionLabelSurface,numThickness,1);

%Get faces for the elements for visualization
[Fp,CFp]=element2patch(Ep,regionLabelElements,'penta6');

Ep_con  = Ep(regionLabelElements==1,:); %Connective tissue
Ep_cell = Ep(regionLabelElements~=1,:); %Cells

%%

cFigure; hold on;
gpatch(Fp,Vp,CFp,'k',0.5);

axisGeom; set(gca,'FontSize',fontSize);
colormap(cMap); icolorbar;
drawnow;

%% Define fibre directions for the solid mesh

V_Ep_con=patchCentre(Ep(regionLabelElements==1,:),Vp);
Np_con1=repmat(N_con1,numThickness,1);

V_Ep_cell=patchCentre(Ep(regionLabelElements~=1,:),Vp);
Np_cell=repmat([0 0 1],size(V_Ep_cell,1),1);

%%

cFigure; hold on;
gpatch(Fp,Vp,CFp,'none',0.05);
quiverLine(V_Ep_con,Np_con1,pointSpacing,'b',lineWidthVec,3,1);
quiverLine(V_Ep_cell,Np_cell,pointSpacing,'r',lineWidthVec,3,1);
colormap(cMap); %icolorbar;
axisGeom; set(gca,'FontSize',fontSize);
colormap(cMap); icolorbar;
drawnow;

%% CREATE FEBio FILE

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

materialName2='Material2';
febio_spec.Material.material{2}.ATTR.name=materialName2;
febio_spec.Material.material{2}.ATTR.type='neo-Hookean';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.E=E_youngs1;
febio_spec.Material.material{2}.v=nu1;

% Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(Vp,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=Vp; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='penta6'; %Element type
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(Ep_cell,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=Ep_cell; %The element matrix

partName2='Part2';
febio_spec.Mesh.Elements{2}.ATTR.name=partName2; %Name of this part
febio_spec.Mesh.Elements{2}.ATTR.type='penta6'; %Element type
febio_spec.Mesh.Elements{2}.elem.ATTR.id=(1:1:size(Ep_con,1))'; %Element id's
febio_spec.Mesh.Elements{2}.elem.VAL=Ep_con; %The element matrix

% % -> NodeSets
% nodeSetName1='bcSupportList_X';
% nodeSetName2='bcSupportList_Y';
% nodeSetName3='bcSupportList_Z';
% nodeSetName4='bcPrescribeList';
% 
% febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
% febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcSupportList_X);
% 
% febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
% febio_spec.Mesh.NodeSet{2}.VAL=mrow(bcSupportList_Y);
% 
% febio_spec.Mesh.NodeSet{3}.ATTR.name=nodeSetName3;
% febio_spec.Mesh.NodeSet{3}.VAL=mrow(bcSupportList_Z);
% 
% febio_spec.Mesh.NodeSet{4}.ATTR.name=nodeSetName4;
% febio_spec.Mesh.NodeSet{4}.VAL=mrow(bcPrescribeList);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain{1}.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat=materialName1;

febio_spec.MeshDomains.SolidDomain{2}.ATTR.name=partName2;
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat=materialName2;

%Boundary condition section 
% -> Fix boundary conditions
% febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_x';
% febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
% febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
% febio_spec.Boundary.bc{1}.x_dof=1;
% febio_spec.Boundary.bc{1}.y_dof=0;
% febio_spec.Boundary.bc{1}.z_dof=0;
% 
% febio_spec.Boundary.bc{2}.ATTR.name='zero_displacement_y';
% febio_spec.Boundary.bc{2}.ATTR.type='zero displacement';
% febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName2;
% febio_spec.Boundary.bc{2}.x_dof=0;
% febio_spec.Boundary.bc{2}.y_dof=1;
% febio_spec.Boundary.bc{2}.z_dof=0;
% 
% febio_spec.Boundary.bc{3}.ATTR.name='zero_displacement_z';
% febio_spec.Boundary.bc{3}.ATTR.type='zero displacement';
% febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName3;
% febio_spec.Boundary.bc{3}.x_dof=0;
% febio_spec.Boundary.bc{3}.y_dof=0;
% febio_spec.Boundary.bc{3}.z_dof=1;
% 
% febio_spec.Boundary.bc{4}.ATTR.name='prescibed_displacement_z';
% febio_spec.Boundary.bc{4}.ATTR.type='prescribed displacement';
% febio_spec.Boundary.bc{4}.ATTR.node_set=nodeSetName4;
% febio_spec.Boundary.bc{4}.dof='z';
% febio_spec.Boundary.bc{4}.value.ATTR.lc=1;
% febio_spec.Boundary.bc{4}.value.VAL=displacementMagnitude;
% febio_spec.Boundary.bc{4}.relative=0;

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

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='sz';
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
%system(['gedit ',febioFebFileName,' &']);


%%


function [E_new,V_new]=resampleEdgeSet(E_now,VC,pointSpacing,smoothWeight)

if smoothWeight>1
    smoothWeight=1;
end

if smoothWeight<0
    smoothWeight=0;
end

if size(E_now,1)>1
    indNow=edgeListToCurve(E_now); %Curve indices for resampling
else
    indNow=E_now;
end

if smoothWeight==0 %Simply linearly upsample cell edges        
    if numel(indNow)==2
            V_new=evenlySpaceCurve(VC(indNow,:),pointSpacing,'linear',0); %Resample
            E_new=[(1:1:size(V_new,1)-1)' (2:1:size(V_new,1))']+size(VC,1);
    else
        if indNow(1)==indNow(end)
            [V_new]=evenlySpaceCurve(VC(indNow(1:end-1),:),pointSpacing,'linear',1,1:1:numel(indNow)-1);
            E_new=[(1:1:size(V_new,1))' [(2:1:size(V_new,1))'; 1]]+size(VC,1);
        else
            [V_new]=evenlySpaceCurve(VC(indNow,:),pointSpacing,'linear',0,1:1:numel(indNow));
            E_new=[(1:1:size(V_new,1)-1)' (2:1:size(V_new,1))']+size(VC,1);
        end
    end

else
    
    if indNow(1)==indNow(end) %Closed curve
        
        V_now=subCurve(VC(indNow(1:end-1),:),1,1);
        
        ind1=(1:2:size(V_now,1))';
        ind2=ind1+1; ind2(ind2>size(V_now,1))=1;
        ind3=ind1-1; ind3(ind3<1)=size(V_now,1);
        
        V_smooth=V_now;
        p=(V_smooth(ind2,:)+V_smooth(ind3,:))/2;
        n=V_smooth(ind1,:)-p;
        c=p-n;
        rr=V_now(ind1,:)-c;
        j=(rr./sqrt(2))+c;%c+nn.*rc;
        V_smooth(ind1,:)=j;
        
        V_smooth=evenlySpaceCurve(V_smooth,pointSpacing/2,'spline',1); %Resample
        V_raw=evenlySampleCurve(V_now,size(V_smooth,1),'linear',1); %Resample        
        V_new=evenlySpaceCurve((1-smoothWeight)*V_raw+(smoothWeight)*V_smooth,pointSpacing,'pchip',1); %Resample
        
        E_new=[(1:1:size(V_new,1))' [(2:1:size(V_new,1))'; 1]]+size(VC,1);
        
    else %Open curve
        
        if numel(indNow)==2
            V_new=evenlySpaceCurve(VC(indNow,:),pointSpacing,'linear',0); %Resample
        else
            
            V_now=subCurve(VC(indNow,:),1,0);
            ind1=(1:2:size(V_now,1))';
            ind1=ind1(ind1~=1 & ind1~=size(V_now,1));
            ind2=ind1+1;
            ind3=ind1-1;
            
            V_smooth=V_now;
            p=(V_smooth(ind2,:)+V_smooth(ind3,:))/2;
            n=V_smooth(ind1,:)-p;
            c=p-n;
            rr=V_now(ind1,:)-c;
            j=(rr./sqrt(2))+c;%c+nn.*rc;
            V_smooth(ind1,:)=j;
            
            V_smooth=evenlySpaceCurve(V_smooth,pointSpacing,'spline',0); %Resample
            V_raw=evenlySampleCurve(V_now,size(V_smooth,1),'linear',0); %Resample            
            V_new=evenlySpaceCurve((1-smoothWeight)*V_raw+(smoothWeight)*V_smooth,pointSpacing,'pchip',0); %Resample
            
        end
        
        E_new=[(1:1:size(V_new,1)-1)' (2:1:size(V_new,1))']+size(VC,1);
        
    end
end

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
