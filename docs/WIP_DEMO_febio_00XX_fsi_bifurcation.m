%% DEMO_febio_0036_sphere_branches_slide_body_force
% Below is a demonstration for:
% 
% * Building geometry for a spherical blob with hexahedral elements
% which is pushed through a branched network. 
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
cMap=blood(250);

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
% febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

pointSpacing=1; 
contourLevel=1;

voxelSize=pointSpacing*ones(1,3); 

optionStructRay.tolEps        = 1e-6;
optionStructRay.triSide       = -1;
optionStructRay.rayType       = 'ray';
optionStructRay.exclusionType = 'inclusive';
optionStructRay.paired        = 0; 

trunkLengthPreSplit=20;
trunkLengthPostSplit=20;
branchLength=20; 
branchAngleDeg=60; 
radiusTrunkInner=6;
radiusBranchInner=4;
wallThickness=1;

%%

radiusTrunkOuter=radiusTrunkInner+wallThickness;
radiusBranchOuter=radiusBranchInner+wallThickness;

pointSpacingCurves=min(voxelSize)/2;

%% Define branch curves

Q=euler2DCM([0 branchAngleDeg/180*pi 0]);
Vc_trunk1 = [0 0 -trunkLengthPreSplit; 0 0 0;];
Vc_trunk1 = evenlySpaceCurve(Vc_trunk1,pointSpacingCurves); 
Vc_trunk2 = [0 0 0; 0 0 trunkLengthPreSplit];
Vc_trunk2 = evenlySpaceCurve(Vc_trunk2,pointSpacingCurves); 
Vc_trunk=[Vc_trunk1;Vc_trunk2(2:end,:)];

Vc_branch = [0 0 0; 0 0 branchLength]*Q;
Vc_branch=evenlySpaceCurve(Vc_branch,pointSpacingCurves); 

Ec_trunk=[(1:1:size(Vc_trunk,1)-1)' (2:1:size(Vc_trunk,1))'];
Rc_trunk=radiusTrunkOuter.*ones(size(Vc_trunk,1),1);

Ec_branch=[(1:1:size(Vc_branch,1)-1)' (2:1:size(Vc_branch,1))'];
Rc_branch=radiusBranchOuter.*ones(size(Vc_trunk,1),1);

R_tree=[Rc_trunk;Rc_branch];
[E_tree,V_tree,C_tree]=joinElementSets({Ec_trunk,Ec_branch},{Vc_trunk,Vc_branch});

[E_tree,V_tree,ind1]=mergeVertices(E_tree,V_tree);
R_tree=R_tree(ind1);

% Q=euler2DCM(rand(1,3)*pi);
% V_tree=V_tree*Q;

%%
%Get branch intersection points
[~,~,~,Ac]=cunique(E_tree);
logicEnd=Ac==1;
indEnd=unique(E_tree(logicEnd));
logicCross=Ac>2;
indCrossing=unique(E_tree(logicCross));

numEdges=size(E_tree,1);
numEdgeVertices=size(E_tree,2);
edgeVertexConnectivity=E_tree;
numVertices=size(V_tree,1);

ind=(1:1:numEdges)';
ind=ind(:,ones(1,numEdgeVertices));
ind=ind(:);
vertexEdgeConnectivity=sparse(edgeVertexConnectivity(:),ind,ind,numVertices,numEdges);
vertexEdgeConnectivity=sort(vertexEdgeConnectivity,2,'descend');
[~,J,~] = find(vertexEdgeConnectivity);
vertexEdgeConnectivity=full(vertexEdgeConnectivity(:,1:max(J)));

c=nan(size(vertexEdgeConnectivity));
c(vertexEdgeConnectivity>0)=C_tree(vertexEdgeConnectivity(vertexEdgeConnectivity>0));

cc=mean(c,2,'omitnan')~=c(:,1);
indCrossingColor=find(cc);


indCrossingColorFalse=indCrossingColor(~ismember(indCrossingColor,indCrossing));
C_branch=C_tree;
for q=1:1:numel(indCrossingColorFalse)
    C_branch(C_branch==c(indCrossingColorFalse(q),2))=c(indCrossingColorFalse(q),1);
end
[~,~,C_branch]=unique(C_branch);

%%
 
cFigure; hold on; 
gpatch(E_tree,V_tree,R_tree,R_tree,0,5);
plotV(V_tree(indEnd,:),'y.','MarkerSize',50);
plotV(V_tree(indCrossing,:),'g.','MarkerSize',50);
axisGeom; camlight headlight;
colormap spectral; colorbar; 
gdrawnow;

%%
rMax=max(R_tree(:));
expandFactor=1.5; 

x=min(V_tree(:,1))-expandFactor*rMax:voxelSize:max(V_tree(:,1))+expandFactor*rMax;
y=min(V_tree(:,2))-expandFactor*rMax:voxelSize:max(V_tree(:,2))+expandFactor*rMax;
z=min(V_tree(:,3))-expandFactor*rMax:voxelSize:max(V_tree(:,3))+expandFactor*rMax;
[X,Y,Z]=ndgrid(x,y,z);
siz=size(X);
V_grid=[X(:) Y(:) Z(:)];

%%
 
cFigure; hold on; 
gpatch(E_tree,V_tree,R_tree,R_tree,0,5);
plotV(V_grid,'k.','markerSize',1);
axisGeom; camlight headlight;
colormap spectral; colorbar; 
gdrawnow;

%%

%Compute all-to-all distances
logicGrid=true(size(V_grid,1),1);
D=distND(V_tree,V_grid(logicGrid,:));

%Compute closest point distances and indices
distType=1;
switch distType
    case 1
        %Compute closest point distances and indices
        [~,I_min]=min(D,[],1); %The index I_min is for rows in D
        J_min=1:1:size(D,2); %Generate column indices for D
        indMin=sub2ind(size(D),I_min,J_min); %Convert to linear indices
        
        %Normalize distances using local radii
        D_norm=D./R_tree(:,ones(nnz(logicGrid),1));
        [D_min_norm,~]=min(D_norm,[],1); %The index I_min is for rows in D
    case 2
        %Normalize distances using local radii
        D_norm=D./R_tree(:,ones(nnz(logicGrid),1));
        [D_min_norm,I_min]=min(D_norm,[],1); %The index I_min is for rows in D
        J_min=1:1:size(D,2); %Generate column indices for D
        indMin=sub2ind(size(D),I_min,J_min); %Convert to linear indices
end

maxDistLevel=2;
M=maxDistLevel.*ones(siz);
M(logicGrid)=D_min_norm;%

%%

vizStruct.colormap=spectral(250);
vizStruct.clim=[0 contourLevel*2];
sv3(M,voxelSize,vizStruct); 
camlight headlight; 

%%

imOrigin=min(V_grid,[],1)-voxelSize/2;

%%

controlPar_isosurface.nSub=[1 1 1];%round(max(v)/2./v);
controlPar_isosurface.capOpt=0; %Option to cap open ended surfaces
controlPar_isosurface.voxelSize=voxelSize;
controlPar_isosurface.contourLevel=contourLevel;
[Fi,Vi]=levelset2isosurface(M,controlPar_isosurface); %Get iso-surface
Fi=fliplr(Fi); %Invert face orientation

Fi_sorted=sort(Fi,2);
logicInvalid=any(diff(Fi_sorted,1,2)==0,2);
Fi=Fi(~logicInvalid,:);
[Fi,Vi]=patchCleanUnused(Fi,Vi);
Vi=Vi(:,[2 1 3]);
Vi=Vi+imOrigin;
[Fi,Vi]=triSurfRemoveThreeConnect(Fi,Vi);
[Fi,Vi]=patchCleanUnused(Fi,Vi);

%% Cut ends
Fi_uncut=Fi;
Vi_uncut=Vi;

for qr=1:1:numel(indEnd)
    
    dCurveEnds=meshDistMarch(E_tree,V_tree,indEnd(qr));
    
    indEdgeMember=find(any(E_tree==indEnd(qr),2));
    c_branch=C_branch(indEdgeMember);
    r_end=R_tree(indEnd(qr));
    
    logicClose=all(dCurveEnds(E_tree)<=r_end,2);
    
    logicEdgeSet=logicClose & C_branch==C_branch(indEdgeMember);
    
    indList=edgeListToCurve(E_tree(logicEdgeSet,:));
    
    if indList(1)~=indEnd(qr)
        indList=indList(end:-1:1); %Invert order
    end
    
    e=E_tree(any(E_tree==indEnd(qr),2),:);    
    v1=V_tree(indEnd(qr),:);
    v2=V_tree(e(e~=indEnd(qr)),:);
    n1=vecnormalize(v1-v2);
        
    V_intersect=triSurfRayTrace(v1,n1,Fi,Vi,optionStructRay);
    
    
    if ~isempty(V_intersect)
        if size(V_intersect,1)>1
            [~,indMin]=minDist(v1,V_intersect);
            V_intersect=V_intersect(indMin,:);
        end
       
        [~,indCloseIntersection]=minDist(V_intersect,Vi);
        
        dFaceEnds=meshDistMarch(Fi,Vi,indCloseIntersection);
        dFaceEnds(isnan(dFaceEnds))=max(dFaceEnds(:));
        
        logicExclude=all(dFaceEnds(Fi)>(r_end*pi),2);
        
        [Fii,Vii,logicExclude_ii,logicSide,Ebi]=triSurfSlice(Fi,Vi,logicExclude,v1,n1,pointSpacing/100,logicExclude);
        Fi=Fii(logicSide | logicExclude_ii,:);
        Vi=Vii;
        [Fi,Vi]=patchCleanUnused(Fi,Vi);
        
%         cFigure; hold on;
%         gpatch(Fi,Vi,'w','k',1);
% %         gpatch(Fi,Vi,logicExclude,'k',1);
% %         gpatch(Fii,Vii+4,logicSide,'r',1);
% %         gpatch(Fii,Vii+8,logicExclude_ii,'g',1);
%         %         patchNormPlot(Fi,Vi)
%         plotV(V_intersect,'r.','MarkerSize',25)
%         quiverVec(v1,n1*r_end,[],'b');
%         plotV(v1,'r.','MarkerSize',50)
%         axisGeom(gca,fontSize); camlight headlight;
%         icolorbar;
%         gdrawnow;
%         
    else
        warning('No intersection found for vessel end cutting');
    end
end

%% Remesh using geomgram
optionStructGG.pointSpacing=pointSpacing;
[Fi,Vi]=ggremesh(Fi,Vi,optionStructGG);

%% Fix boundary curve coordinates
% The ggremesh function has a smoothing effect and associated boundary
% "shrinkage". This code pushes the boundaries back. 

Eb=patchBoundary(Fi); %All boundary edges

%Use grouping to seperate
optionStructGroup.outputType='label';
G=tesgroup(Eb,optionStructGroup);

%Loop over all groups and correct
for q=1:1:max(G(:))
    eNow=Eb(G==q,:); %The current edge set 
    indNow=unique(eNow); %The current vertex indices
    Vm=mean(Vi(indNow,:),1); %The mean or centre of the current boundary curve
    [~,indMin]=minDist(Vm,V_tree(indEnd,:)); %Index for nearest centre line end point
    indTreeNow=indEnd(indMin); %Current index for tree end point
    logicTreeEdgeNow=any(E_tree==indTreeNow,2); %Logic for edge on tree attached to end point
    
    e=E_tree(logicTreeEdgeNow,:); %Get edge at end point
    v1=V_tree(indTreeNow,:); %Get end point 
    v2=V_tree(e(e~=indTreeNow),:); %Get "previous" point to form direction vector
    n1=vecnormalize(v1-v2); %Tree end direction vectors    
    R=vecPair2Rot(n1,[0 0 1]); %Contruct rotation matrix
    
    %Rotate 
    Vi=Vi-v1; %Centre on end point
    Vi=Vi*R'; %Rotate to XY plane
    Vi(indNow,3)=0; %Force Z back to zero after rotation/centering
    Vi=Vi*R; %Rotate back
    Vi=Vi+v1; %Translate back
end

%%
% Visualize mesh 

cFigure; hold on; 
gpatch(Fi,Vi,'w','k');
% gpatch(Fi_uncut,Vi_uncut,'bw','none',0.5);
% patchNormPlot(Fi,Vi);
gpatch(E_tree,V_tree,R_tree,R_tree,0,5);
plotV(V_tree(indEnd,:),'y.','MarkerSize',50);
plotV(V_tree(indCrossing,:),'g.','MarkerSize',50);
axisGeom; camlight headlight;
colormap spectral; colorbar; 
gdrawnow;

%% Thicken into pentrahedral elements

numSteps=2;
dirSet=-1;
[Ep,Vp,Fp1,Fp2]=patchThick(Fi,Vi,dirSet,wallThickness,numSteps);

Fp=element2patch(Ep,[],'penta6');

%%
% Visualize mesh 

cFigure; hold on; 
gpatch(Fp,Vp,'w','k',0.5);
gpatch(Fp1,Vp,'rw','k');
gpatch(Fp2,Vp,'bw','k');
axisGeom; camlight headlight;
gdrawnow;

%% Cap all ends

Eb=patchBoundary(Fp2); %All boundary edges

%Use grouping to seperate
optionStructGroup.outputType='label';
G=tesgroup(Eb,optionStructGroup);

%Loop over all groups and correct
FT=cell(max(G(:)),1);
VT=cell(max(G(:)),1);

for q=1:1:max(G(:))
    
    eNow=Eb(G==q,:); %The current edge set 
    indNow=edgeListToCurve(eNow); %The current vertex indices
    indNow=indNow(1:end-1); 
    Vm=mean(Vp(indNow,:),1); %The mean or centre of the current boundary curve
    [~,indMin]=minDist(Vm,V_tree(indEnd,:)); %Index for nearest centre line end point
    indTreeNow=indEnd(indMin); %Current index for tree end point
    logicTreeEdgeNow=any(E_tree==indTreeNow,2); %Logic for edge on tree attached to end point
    
    e=E_tree(logicTreeEdgeNow,:); %Get edge at end point
    v1=V_tree(indTreeNow,:); %Get end point 
    v2=V_tree(e(e~=indTreeNow),:); %Get "previous" point to form direction vector
    n1=vecnormalize(v1-v2); %Tree end direction vectors    
    R=vecPair2Rot(n1,[0 0 1]); %Contruct rotation matrix
    
    
    Vc=(Vp(indNow,:)-v1)*R';
    
    inputStructureRegionTriMesh.regionCell={Vc(:,[1 2])};
    inputStructureRegionTriMesh.pointSpacing=pointSpacing;
    inputStructureRegionTriMesh.resampleCurveOpt=0; %Turn on/off curve resampling
    inputStructureRegionTriMesh.plotOn=0;
    [Ft,Vt]=regionTriMesh2D(inputStructureRegionTriMesh);    
    Vt(:,3)=0;
    FT{q}=Ft;
    VT{q}=(Vt*R)+v1;
    
end

%%

cFigure; hold on; 
gpatch(FT,VT,'y','k',1);
patchNormPlot(FT,VT);
gpatch(Fp,Vp,'w','k',0.5);
gpatch(Fp1,Vp,'rw','k');
gpatch(Fp2,Vp,'bw','k');
axisGeom; camlight headlight;
gdrawnow;


%%

%%

cFigure; hold on; 
title('Fluid domain')
gpatch(FT,VT,'rw','k',1);
patchNormPlot(FT,VT);
gpatch(Fp2,Vp,'bw','k');
patchNormPlot(Fp2,Vp);

axisGeom; camlight headlight;
gdrawnow;


fdasdff


%%

%Control settings
cPar.sphereRadius=sphereRadius;
cPar.coreRadius=sphereRadius.*0.75;
cPar.numElementsCore=ceil(sphereRadius/pointSpacing); 
cPar.numElementsMantel=ceil((sphereRadius-cPar.coreRadius)/pointSpacing); 
cPar.makeHollow=0;
cPar.outputStructType=2;

%Creating sphere
[meshStruct]=hexMeshSphere(cPar);

%Access ouput
E_blob=meshStruct.elements; %The elements 
V_blob=meshStruct.nodes; %The vertices
Fb_blob=meshStruct.facesBoundary; %The boundary faces

%%
%Create cut view

hFig=cFigure;
subplot(1,2,1); hold on;
title('The hexahedral mesh sphere','FontSize',fontSize);
gpatch(Fb_blob,V_blob,'r');
axisGeom(gca,fontSize);
camlight headlight;

hs=subplot(1,2,2); hold on; 
title('Cut view of solid mesh','FontSize',fontSize);
optionStruct.hFig=[hFig hs];
gpatch(Fb_blob,V_blob,'kw','none',0.25);
meshView(meshStruct,optionStruct);
axisGeom(gca,fontSize);
gdrawnow;

gdrawnow;

%% Shift branch set

Vi=Vi(:,[2 1 3]);
Vi=Vi+imOrigin;

%% Join model node sets

V=[V_blob; Vi; ];
F_tube=Fi+size(V_blob,1);
F_tube=fliplr(F_tube);

center_of_mass_tube=mean(Vi,1);

%%
% Visualizing model

cFigure; hold on;
gtitle('Model components',fontSize);
hl(1)=gpatch(Fb_blob,V,'rw','k',0.8);
hl(2)=gpatch(F_tube,V,'bw','k',0.8);
patchNormPlot(F_tube,V);
legend(hl,{'Clot','Vessel'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 
gdrawnow;

%% Get contact surfaces
%

F_contact_blob=Fb_blob;

%%
% Visualize contact surfaces

cFigure; hold on;
title('Tube blob contact pair','fontsize',fontSize);
hl(1)=gpatch(F_tube,V,'rw','k',1);
patchNormPlot(F_tube,V);
hl(2)=gpatch(F_contact_blob,V,'gw','k',1);
patchNormPlot(F_contact_blob,V);
legend(hl,{'Secondary','Primary'}); clear hl;
axisGeom(gca,fontSize);
camlight headlight; 
gdrawnow;

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
switch materialType
    case 0
        febio_spec.Material.material{1}.ATTR.type='Ogden unconstrained';
        febio_spec.Material.material{1}.ATTR.id=1;
        febio_spec.Material.material{1}.c1=c1;
        febio_spec.Material.material{1}.m1=m1;
        febio_spec.Material.material{1}.c2=c1;
        febio_spec.Material.material{1}.m2=-m1;
        febio_spec.Material.material{1}.cp=k;
        febio_spec.Material.material{1}.density=materialDensity;
    case 1
        %Viscoelastic part
        febio_spec.Material.material{1}.ATTR.type='viscoelastic';        
        febio_spec.Material.material{1}.ATTR.id=1;
        febio_spec.Material.material{1}.g1=g1;
        febio_spec.Material.material{1}.t1=t1;
        febio_spec.Material.material{1}.density=materialDensity;
        
        %Elastic part
        febio_spec.Material.material{1}.elastic{1}.ATTR.type='Ogden unconstrained';
        febio_spec.Material.material{1}.elastic{1}.c1=c1;
        febio_spec.Material.material{1}.elastic{1}.m1=m1;
        febio_spec.Material.material{1}.elastic{1}.c2=c1;
        febio_spec.Material.material{1}.elastic{1}.m2=-m1;
        febio_spec.Material.material{1}.elastic{1}.cp=k;
        febio_spec.Material.material{1}.elastic{1}.density=materialDensity;
end

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
surfaceName1='contactSurface1';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(F_tube,1))';
febio_spec.Mesh.Surface{1}.quad4.VAL=F_tube;

surfaceName2='contactSurface2';
febio_spec.Mesh.Surface{2}.ATTR.name=surfaceName2;
febio_spec.Mesh.Surface{2}.quad4.ATTR.id=(1:1:size(F_contact_blob,1))';
febio_spec.Mesh.Surface{2}.quad4.VAL=F_contact_blob;

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

%Loads section
% -> Body load
febio_spec.Loads.body_load.ATTR.type='const';
febio_spec.Loads.body_load.x.ATTR.lc=1;
febio_spec.Loads.body_load.x.VAL=bodyLoadMagnitude;
      

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
febio_spec.Contact.contact{1}.search_radius=0.1*sqrt(sum((max(V,[],1)-min(V,[],1)).^2,2)); 
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
febioAnalysis.runMode='internal';%'internal';

[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%% Import FEBio results 

if 1%runFlag==1 %i.e. a succesful run
    
       
    %%
    % Importing nodal displacements from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
    
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
    hf=cFigure; hold on;
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp=gpatch(Fb_blob,V_DEF(:,:,end),DN_magnitude,'k',1); %Add graphics object to animate
    hp.FaceColor='interp';
    gpatch(F_tube,V,'w','none',0.5); %Add graphics object to animate
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]/3);
    axis(axisLim(V_DEF)); %Set axis limits statically  
    camlight headlight;
    gdrawnow;
    
    % Set up animation features
    animStruct.Time=timeVec; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments
        DN_magnitude=sqrt(sum(N_disp_mat(:,:,qt).^2,2)); %Current displacement magnitue
        
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
